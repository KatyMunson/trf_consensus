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

Nothing is dropped — every original cluster call still appears in the
final outputs, but non-representative members of a family are marked
is_representative=False and demoted in the RepeatMasker library's
classification field, rather than removed.

Outputs:
  --out-lib-fasta: EVERY cluster's consensus, headers as
      >{final_name}__{entry_id}#{classification} for the representative, and
      >{final_name}__{entry_id}#{classification}/redundant_with_{final_name}
      for other family members. Sorted family-discovery-order, representative
      first within each family.
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
        "final_name", "motif_length", "is_representative", "source_entry_id",
        "source_cluster_label", "source_individual", "source_assembly_method",
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

    ranked = sorted(clusters, key=lambda c: float(c["total_copy_number"]), reverse=True)

    winners = []                       # winner cluster_uids, in discovery order
    redundant_with = {}                # loser cluster_uid -> winner cluster_uid
    group_members = defaultdict(list)  # winner cluster_uid -> [winner, loser, loser, ...]

    for c in ranked:
        uid = c["cluster_uid"]
        match = None
        for w in winners:
            if w in flagged[uid]:
                match = w
                break  # winners is already support-sorted; first hit is the strongest
        if match is None:
            winners.append(uid)
            group_members[uid].append(uid)
        else:
            redundant_with[uid] = match
            group_members[match].append(uid)

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

    # --- final summary table, family-discovery order, representative first within each family ---
    rows_out = []
    for w in winners:
        final_name = final_name_of_winner[w]
        motif_length = int(by_uid[w]["consensus_length"])
        prov = provenance[w]
        for uid in group_members[w]:
            c = by_uid[uid]
            rows_out.append({
                "final_name": final_name, "motif_length": motif_length,
                "is_representative": uid == w,
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
    # worked example never exercises.
    seen_header_keys = set()
    for row in rows_out:
        key = (row["final_name"], row["source_entry_id"])
        if key in seen_header_keys:
            sys.exit(
                f"ERROR: family '{row['final_name']}' has more than one cluster from entry "
                f"'{row['source_entry_id']}' — the {{final_name}}__{{entry_id}} library header "
                f"format can't disambiguate them. Resolve manually (this is not expected to "
                f"happen for one cluster-per-bin-per-entry input)."
            )
        seen_header_keys.add(key)

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
        f"{len(winners)} families, {n_redundant} non-representative members\n"
    )


if __name__ == "__main__":
    main()
