#!/usr/bin/env python
"""
resolve_redundancy.py

Takes the master cluster table (all_clusters.tsv), the pairwise comparison
(cross_cluster_comparison.tsv: flag_similar / flag_multiple_of), and the
combined consensus FASTA, and resolves redundancy.

IMPORTANT — this does NOT use transitive union-find over flagged pairs.
An earlier version did, and it produced meaningless mega-groups: if A~B is
flagged and B~C is flagged, plain union-find merges A, B, and C into one
group even when A and C have no direct relationship at all — and chaining
this way across a few hundred clusters can (and did, on real data) collapse
dozens of genuinely unrelated candidates spanning an enormous period range
into one "family," with no cluster in the group ever compared to most of
the others.

Instead: process clusters in descending support order (highest
n_input_sequences first). Each cluster is checked ONLY against clusters
already confirmed as group winners so far — never against another loser,
never chained through an intermediate. If it's flagged as similar to or a
multiple of an existing winner, it's marked redundant with THAT winner
specifically (a direct, individually-verified relationship). Otherwise it
becomes a new winner itself. This is the same non-chained pattern used to
fix the equivalent bug in the old triage_candidates.py.

Nothing is dropped — every cluster still appears in the final outputs, but
losers are marked redundant and demoted (both in the summary table and in
the RepeatMasker library's classification field) rather than removed.

Outputs:
  --out-lib-fasta: EVERY cluster's consensus, headers as
      >{label}#{classification} for winners, and
      >{label}#{classification}/redundant_with_{winner} for demoted ones.
      Sorted winners-first, then by support.
  --out-summary: all_clusters.tsv plus is_redundant / redundant_with /
      group_size columns — the final master QC table.
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


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--clusters-tsv", required=True, help="all_clusters.tsv")
    ap.add_argument("--comparison-tsv", required=True, help="cross_cluster_comparison.tsv")
    ap.add_argument("--consensus-fasta", required=True, help="all_clusters_consensus.fasta")
    ap.add_argument("--classification", default="Satellite")
    ap.add_argument("--out-lib-fasta", required=True)
    ap.add_argument("--out-summary", required=True)
    args = ap.parse_args()

    clusters = read_tsv(args.clusters_tsv)
    comparisons = read_tsv(args.comparison_tsv)
    seqs = read_multi_fasta(args.consensus_fasta)

    by_label = {c["label"]: c for c in clusters}

    # flagged[label] -> set of OTHER labels it's flagged similar/multiple-of
    flagged = defaultdict(set)
    n_edges = 0
    for r in comparisons:
        if r.get("flag_similar") == "True" or r.get("flag_multiple_of") == "True":
            a, b = r["label_A"], r["label_B"]
            if a in by_label and b in by_label:
                flagged[a].add(b)
                flagged[b].add(a)
                n_edges += 1

    ranked = sorted(clusters, key=lambda c: int(c["n_input_sequences"]), reverse=True)

    winners = []          # labels confirmed as winners so far, in order found
    redundant_with = {}   # label -> winner label
    group_members = defaultdict(list)  # winner label -> list of labels redundant with it (incl. itself)

    for c in ranked:
        label = c["label"]
        match = None
        for w in winners:
            if w in flagged[label]:
                match = w
                break  # winners list is already support-sorted; first hit is the strongest
        if match is None:
            winners.append(label)
            group_members[label].append(label)
        else:
            redundant_with[label] = match
            group_members[match].append(label)

    group_size_of = {}
    for winner, members in group_members.items():
        for m in members:
            group_size_of[m] = len(members)

    # --- final summary table ---
    out_columns = ["label", "motif", "target_period", "window", "rank",
                   "n_input_sequences", "pct_of_bin", "consensus_length", "gc_content",
                   "is_redundant", "redundant_with", "group_size"]
    rows_out = []
    for c in clusters:
        label = c["label"]
        is_redundant = label in redundant_with
        row = dict(c)
        row["is_redundant"] = is_redundant
        row["redundant_with"] = redundant_with.get(label, "NA")
        row["group_size"] = group_size_of[label]
        rows_out.append(row)

    # Winners first, then by support descending
    rows_out.sort(key=lambda r: (r["is_redundant"], -int(r["n_input_sequences"])))

    with open(args.out_summary, "w") as out:
        out.write("\t".join(out_columns) + "\n")
        for row in rows_out:
            out.write("\t".join(str(row[c]) for c in out_columns) + "\n")

    # --- final RepeatMasker library ---
    with open(args.out_lib_fasta, "w") as out:
        for row in rows_out:
            label = row["label"]
            seq = seqs.get(label, "")
            if not seq:
                continue
            if row["is_redundant"]:
                classification = f"{args.classification}/redundant_with_{row['redundant_with']}"
            else:
                classification = args.classification
            out.write(f">{label}#{classification}\n{seq}\n")

    n_redundant = sum(1 for r in rows_out if r["is_redundant"])
    max_group = max(group_size_of.values()) if group_size_of else 0
    sys.stderr.write(
        f"[resolve_redundancy] {len(clusters)} clusters, {n_edges} flagged relationships, "
        f"{len(winners)} winners, {n_redundant} marked redundant, largest group={max_group}\n"
    )


if __name__ == "__main__":
    main()
