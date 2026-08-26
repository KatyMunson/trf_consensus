#!/usr/bin/env python
"""
resolve_redundancy.py

Takes the pooled cross-entry cluster table (filtered_pooled_clusters.tsv),
the pairwise comparison (cross_cluster_comparison.tsv: flag_similar /
flag_multiple_of, keyed on cluster_uid now that clusters are pooled across
entries), and the pooled consensus FASTA, and resolves redundancy into
named families.

IMPORTANT — this does NOT use transitive union-find over flagged pairs.
An earlier version did, and it produced meaningless mega-groups: if A~B is
flagged and B~C is flagged, plain union-find merges A, B, and C into one
group even when A and C have no direct relationship at all — and chaining
this way across a few hundred clusters can (and did, on real data) collapse
dozens of genuinely unrelated candidates spanning an enormous period range
into one "family," with no cluster in the group ever compared to most of
the others.

Instead: process clusters in descending total_copy_number order (the
per-entry summed TRF copy number backing each cluster, not raw locus
count — see rank_family_clusters.py). Each cluster is checked ONLY against
clusters already confirmed as group winners so far — never against
another loser, never chained through an intermediate. If it's flagged as
similar to or a multiple of an existing winner, it joins that winner's
family; otherwise it becomes a new winner (and the representative of a new
family) itself. Same non-chained pattern as before.

TWO-LEVEL RESOLUTION — a within-entry consolidation pre-pass runs before
the global pass described above, using the exact same non-chained
algorithm (see resolve_winners()), just scoped to one entry's own clusters
at a time. This exists because one entry's TRF period scan can nominate
several candidate bins for the same real monomer (e.g. periods 411/417/423
all mutually flagged similar — jitter that survived nms_radius
suppression), and the global pass alone only ever compares each cluster
against already-confirmed GLOBAL winners, never against another candidate
from its own entry. Two such same-entry clusters could therefore each
independently get flagged against the same external winner and land in
the same family without ever being compared to each other — producing two
rows that collide on the {final_name}__{entry_id} library header identity
downstream. The pre-pass fixes this at the source: within each entry,
clusters are ranked and resolved into "sub-winners" first (one sub-winner
per real family present in that entry); only sub-winners are fed into the
unchanged global pass. Every cluster consolidated away within its entry is
still kept in the final outputs (is_representative=False, same as any
other non-representative row) and records which sub-winner it was folded
into via within_entry_consolidated_into (NA for sub-winners themselves).
A sub-winner's total_copy_number is used as-is for global-pass ranking —
deliberately NOT summed across what it consolidated, since overlapping
bin extraction windows (period +/- window) mean the same raw TRF loci
could plausibly have been double-counted into more than one bin already.

The pre-pass alone is not sufficient: two same-entry sub-winners can also
converge on the same family INDIRECTLY, via independent links to a shared
EXTERNAL winner, without ever being directly flagged against each other —
confirmed on real data (two clusters each >= threshold against a common
external winner, but below min_coverage against each other, so
flag_similar/flag_multiple_of both False for the direct pair). The
pre-pass has no visibility into this: it only ever compares candidates
within one entry, never against a shared external target that's only
discovered once the global pass runs. So a second, POST-GLOBAL
deduplication pass runs immediately after the global resolve_winners()
call: within each winner's group_members[w], if more than one member
shares an entry_id (only possible via this indirect-convergence path,
since the pre-pass already ruled out direct same-entry duplicates), the
highest-total_copy_number one is kept and the rest are re-homed under it
(along with anything they'd already consolidated within their own entry).
`w` itself is always correctly kept when it's the duplicated entry, since
resolve_winners()'s ranked processing guarantees w has the highest copy
number among everything in its own group. Together, the two passes give a
structurally exhaustive guarantee: group_members[w] is the complete
definition of family membership, and there are exactly two ways a cluster
joins it (pre-pass within-entry consolidation, or direct/global match
against w) — both are now deduplicated by entry_id, so no two same-entry
clusters can reach the final output claiming the same family, direct or
indirect.

Naming: every newly-confirmed winner gets final_name = "SAT{motif_length}_
{letter}", where motif_length is the WINNER'S OWN exact consensus_length
(no rounding/averaging — different entries can legitimately disagree by a
base or two, since each ran its own independent MSA/consensus over a
different input pool). The letter disambiguates families that happen to
share a motif_length: assigned in winner-discovery order (i.e. descending
total_copy_number among winners), per-length, starting at 'a'.

Provenance columns (n_entries_found/entries_found, n_individuals_present/
individuals_found, n_methods_confirming/methods_found) are computed once
per family and repeated on every row belonging to it. entries_found is
sorted(set(entry_id)); individuals_found/methods_found are derived by
walking entries_found in that same sorted order and deduping each entry's
individual/method while preserving first-occurrence order (not an
alphabetical sort of the values themselves).

Nothing is dropped from --out-summary — every original cluster call still
appears there, but non-representative members of a family are marked
is_representative=False, rather than removed. within_entry_consolidated_into
gives that demotion a second, finer-grained provenance dimension:
is_representative/final_name describe global, cross-entry demotion, while
within_entry_consolidated_into records the (separate, earlier)
within-entry fold, so a row can be within_entry_consolidated_into a
sub-winner that itself later becomes globally is_representative=False.

--out-lib-fasta is narrower: it only ever gets one sequence per entry per
family (the entry's own sub-winner, win or lose globally), since that's
all the {final_name}__{entry_id} header format can uniquely address and a
within-entry-consolidated cluster's sequence is, by construction, already
redundant with its own entry's sub-winner. --out-summary is the complete
record of where every original candidate ended up.

Outputs:
  --out-lib-fasta: one sequence per entry per family (each entry's
      sub-winner), headers as
      >{final_name}__{entry_id}#{classification} for the global representative, and
      >{final_name}__{entry_id}#{classification}/redundant_with_{final_name}
      for every other entry's sub-winner in that family. Sorted
      family-discovery-order, representative first within each family.
  --out-summary: one row per original cluster call — the final harmonized
      master QC table.
"""
import argparse
import sys
from collections import defaultdict


def read_tsv(path):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = [dict(zip(header, line.rstrip("\n").split("\t"))) for line in f if line.strip()]
    return rows


def read_multi_fasta(path):
    seqs = {}
    header, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(seq)
                header = line[1:].split()[0] if line[1:].split() else line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if header is not None:
            seqs[header] = "".join(seq)
    return seqs


def letter_for_index(n):
    """Bijective base-26 letters for n=0,1,2,...: a, b, ..., z, aa, ab, ...,
    az, ba, ... — so a motif_length shared by more than 26 families doesn't
    IndexError, however unlikely that is for real satellite data."""
    n += 1
    s = ""
    while n > 0:
        n -= 1
        s = chr(ord("a") + n % 26) + s
        n //= 26
    return s


def resolve_winners(candidate_uids, flagged, by_uid):
    """Rank candidate_uids by total_copy_number descending; each is checked
    only against already-confirmed winners from THIS SAME candidate_uids
    pool — never against a loser, never chained through an intermediate,
    and never against a cluster outside candidate_uids (this is what scopes
    one call to a single entry vs. globally: `winners` is only ever
    populated from what's passed in, so a match can never reach outside the
    pool). Returns (winners, group_members, redundant_with):
      winners        -- winner cluster_uids, in discovery order
      group_members  -- winner cluster_uid -> [winner, loser, loser, ...]
      redundant_with -- loser cluster_uid -> the winner it joined
    """
    ranked = sorted(candidate_uids, key=lambda u: float(by_uid[u]["total_copy_number"]), reverse=True)
    winners = []
    redundant_with = {}
    group_members = defaultdict(list)
    for uid in ranked:
        match = next((w for w in winners if w in flagged[uid]), None)
        if match is None:
            winners.append(uid)
            group_members[uid].append(uid)
        else:
            redundant_with[uid] = match
            group_members[match].append(uid)
    return winners, group_members, redundant_with


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--clusters-tsv", required=True, help="filtered_pooled_clusters.tsv")
    ap.add_argument("--comparison-tsv", required=True, help="cross_cluster_comparison.tsv")
    ap.add_argument("--consensus-fasta", required=True, help="pooled_consensus.fasta")
    ap.add_argument("--classification", default="Satellite")
    ap.add_argument("--out-lib-fasta", required=True)
    ap.add_argument("--out-summary", required=True)
    args = ap.parse_args()

    clusters = read_tsv(args.clusters_tsv)
    comparisons = read_tsv(args.comparison_tsv)
    seqs = read_multi_fasta(args.consensus_fasta)

    out_columns = [
        "final_name", "motif_length", "is_representative", "within_entry_consolidated_into",
        "source_entry_id", "source_cluster_label", "source_individual", "source_assembly_method",
        "source_phasing_status", "source_copy_number", "source_consensus_length",
        "source_gc_content", "n_entries_found", "entries_found", "n_individuals_present",
        "individuals_found", "n_methods_confirming", "methods_found",
    ]

    if not clusters:
        sys.stderr.write(
            "[resolve_redundancy] no clusters in --clusters-tsv (everything filtered out "
            "upstream?) — writing header-only outputs\n"
        )
        with open(args.out_summary, "w") as out:
            out.write("\t".join(out_columns) + "\n")
        open(args.out_lib_fasta, "w").close()
        return

    by_uid = {c["cluster_uid"]: c for c in clusters}

    # flagged[uid] -> set of OTHER cluster_uids it's flagged similar/multiple-of
    flagged = defaultdict(set)
    n_edges = 0
    for r in comparisons:
        if r.get("flag_similar") == "True" or r.get("flag_multiple_of") == "True":
            a, b = r["label_A"], r["label_B"]
            if a in by_uid and b in by_uid:
                flagged[a].add(b)
                flagged[b].add(a)
                n_edges += 1

    # --- within-entry consolidation pre-pass: same non-chained algorithm,
    # scoped to one entry's own clusters at a time, so two same-entry
    # clusters that both independently match the same external winner get
    # compared to EACH OTHER first, before either ever reaches the global
    # pass (see module docstring's "TWO-LEVEL RESOLUTION" section). ---
    all_subwinners = []                     # every entry's sub-winners, flattened
    entry_group_members = {}                # sub-winner uid -> [sub-winner, consolidated, ...]
    entry_ids = sorted({c["entry_id"] for c in clusters})
    for entry_id in entry_ids:
        entry_uids = [c["cluster_uid"] for c in clusters if c["entry_id"] == entry_id]
        sub_winners, sub_group_members, _ = resolve_winners(entry_uids, flagged, by_uid)
        all_subwinners.extend(sub_winners)
        entry_group_members.update(sub_group_members)

    # --- global pass: unchanged algorithm, now run over sub-winners only.
    # Each sub-winner's total_copy_number is used as-is (not summed across
    # what it consolidated) -- see module docstring. ---
    winners, group_members, redundant_with = resolve_winners(all_subwinners, flagged, by_uid)

    # --- post-global consolidation: within any single family, at most one
    # sub-winner may represent a given entry (the {final_name}__{entry_id}
    # header format requires this). Two same-entry sub-winners can converge
    # on the same family via independent links to a shared external winner
    # without ever being directly flagged against each other -- the
    # within-entry pre-pass above cannot catch this by design (it only ever
    # compares candidates within one entry, never against a shared external
    # target discovered later in this global pass). Resolved here as a
    # final tie-break: for every family, if more than one member sub-winner
    # shares an entry_id, keep the highest-total_copy_number one and
    # re-home the rest -- and everything THEY had already consolidated
    # within their own entry -- under it. `w` is always the keeper when
    # it's the duplicated entry, since resolve_winners() guarantees w has
    # the highest copy number in its own group (every other member was
    # necessarily processed, and so ranked lower, after w was confirmed). ---
    for w in winners:
        by_entry = defaultdict(list)
        for uid in group_members[w]:
            by_entry[by_uid[uid]["entry_id"]].append(uid)
        for entry_id, uids in by_entry.items():
            if len(uids) <= 1:
                continue
            uids_ranked = sorted(uids, key=lambda u: float(by_uid[u]["total_copy_number"]), reverse=True)
            keep, extras = uids_ranked[0], uids_ranked[1:]
            for extra in extras:
                group_members[w].remove(extra)
                entry_group_members[keep].extend(entry_group_members.pop(extra))

    # Computed after both consolidation passes so a cluster only folded in
    # by the post-global pass above (previously its own top-level
    # sub-winner) is correctly counted as consolidated within-entry.
    n_within_entry_consolidated = sum(len(members) - 1 for members in entry_group_members.values())

    # --- naming: one final_name per winner, letters disambiguate shared motif_length ---
    letter_counter = defaultdict(int)
    final_name_of_winner = {}
    for w in winners:
        motif_length = int(by_uid[w]["consensus_length"])
        letter = letter_for_index(letter_counter[motif_length])
        letter_counter[motif_length] += 1
        final_name_of_winner[w] = f"SAT{motif_length}_{letter}"

    # --- per-family provenance, computed once per winner ---
    provenance = {}
    for w in winners:
        members = group_members[w]
        entries_found = sorted({by_uid[m]["entry_id"] for m in members})

        def dedup_in_entry_order(field):
            seen = []
            for entry_id in entries_found:
                value = next(by_uid[m][field] for m in members if by_uid[m]["entry_id"] == entry_id)
                if value not in seen:
                    seen.append(value)
            return seen

        individuals_found = dedup_in_entry_order("individual")
        methods_found = dedup_in_entry_order("assembly_method")
        provenance[w] = {
            "n_entries_found": len(entries_found), "entries_found": ";".join(entries_found),
            "n_individuals_present": len(individuals_found),
            "individuals_found": ";".join(individuals_found),
            "n_methods_confirming": len(methods_found), "methods_found": ";".join(methods_found),
        }

    # --- final summary table, family-discovery order, representative first
    # within each family. Each winner's group_members[w] holds sub-winner
    # uids (at most one per contributing entry); expand each sub-winner
    # back out to every cluster within-entry-consolidated into it, so every
    # original cluster call still gets exactly one row (nothing dropped). ---
    rows_out = []
    for w in winners:
        final_name = final_name_of_winner[w]
        motif_length = int(by_uid[w]["consensus_length"])
        prov = provenance[w]
        for sub_uid in group_members[w]:
            for uid in entry_group_members[sub_uid]:
                c = by_uid[uid]
                within_entry_consolidated_into = (
                    "NA" if uid == sub_uid else by_uid[sub_uid]["source_cluster_label"]
                )
                rows_out.append({
                    "final_name": final_name, "motif_length": motif_length,
                    "is_representative": uid == w,
                    "within_entry_consolidated_into": within_entry_consolidated_into,
                    "source_entry_id": c["entry_id"], "source_cluster_label": c["source_cluster_label"],
                    "source_individual": c["individual"], "source_assembly_method": c["assembly_method"],
                    "source_phasing_status": c["phasing_status"],
                    "source_copy_number": c["total_copy_number"],
                    "source_consensus_length": c["consensus_length"],
                    "source_gc_content": c["gc_content"],
                    **prov,
                })

    # A family with two rows from the same entry_id would collide on the
    # library FASTA's header ({final_name}__{entry_id}); fail loudly rather
    # than invent an unverified disambiguation scheme for a case the spec's
    # worked example never exercises. Scoped to within_entry_consolidated_into
    # == "NA" (i.e. entry sub-winners, the only rows that actually produce a
    # FASTA header) -- a within-entry-consolidated row deliberately shares
    # (final_name, source_entry_id) with its own sub-winner, that's expected
    # and not a collision. After the within-entry consolidation pre-pass
    # this is expected to be unreachable in normal operation (two same-entry
    # sub-winners can only both reach here if they were never flagged
    # against each other) -- if it fires, that's a signal worth investigating
    # on its own, not evidence this fix is incomplete.
    seen_header_keys = {}
    for row in rows_out:
        if row["within_entry_consolidated_into"] != "NA":
            continue
        key = (row["final_name"], row["source_entry_id"])
        if key in seen_header_keys:
            prev = seen_header_keys[key]
            sys.exit(
                f"ERROR: family '{row['final_name']}' has more than one cluster from entry "
                f"'{row['source_entry_id']}':\n"
                f"  {prev['source_cluster_label']} (consensus_length={prev['source_consensus_length']}, "
                f"copy_number={prev['source_copy_number']})\n"
                f"  {row['source_cluster_label']} (consensus_length={row['source_consensus_length']}, "
                f"copy_number={row['source_copy_number']})\n"
                f"Resolve manually."
            )
        seen_header_keys[key] = row



    with open(args.out_summary, "w") as out:
        out.write("\t".join(out_columns) + "\n")
        for row in rows_out:
            out.write("\t".join(str(row[c]) for c in out_columns) + "\n")

    # --- final RepeatMasker library ---
    with open(args.out_lib_fasta, "w") as out:
        for w in winners:
            final_name = final_name_of_winner[w]
            for uid in group_members[w]:
                c = by_uid[uid]
                seq = seqs.get(uid, "")
                if not seq:
                    continue
                header_name = f"{final_name}__{c['entry_id']}"
                if uid == w:
                    classification = args.classification
                else:
                    classification = f"{args.classification}/redundant_with_{final_name}"
                out.write(f">{header_name}#{classification}\n{seq}\n")

    n_redundant = sum(1 for r in rows_out if not r["is_representative"])
    sys.stderr.write(
        f"[resolve_redundancy] {len(clusters)} clusters, {n_edges} flagged relationships, "
        f"{n_within_entry_consolidated} consolidated within-entry ({len(all_subwinners)} entry "
        f"sub-winners fed to the global pass), {len(winners)} families, "
        f"{n_redundant} non-representative members\n"
    )


if __name__ == "__main__":
    main()
