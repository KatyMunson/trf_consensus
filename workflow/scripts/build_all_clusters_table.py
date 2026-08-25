#!/usr/bin/env python
"""
build_all_clusters_table.py

Aggregate every motif's ranked_clusters_summary.tsv (from
rank_family_clusters.py) into one master table, one row per cluster across
every bin scanned — adding the bin's target_period/window for context.
"""
import argparse
import sys


def read_tsv(path):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = [dict(zip(header, line.rstrip("\n").split("\t"))) for line in f if line.strip()]
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--motifs", nargs="+", required=True)
    ap.add_argument("--periods", nargs="+", type=int, required=True)
    ap.add_argument("--windows", nargs="+", type=int, required=True)
    ap.add_argument("--summaries", nargs="+", required=True, help="ranked_clusters_summary.tsv per motif")
    ap.add_argument("--out-tsv", required=True)
    args = ap.parse_args()

    n = len(args.motifs)
    for name, lst in [("periods", args.periods), ("windows", args.windows), ("summaries", args.summaries)]:
        if len(lst) != n:
            sys.exit(f"ERROR: --{name} has {len(lst)} entries, expected {n} (one per motif)")

    columns = ["label", "motif", "target_period", "window", "rank",
               "n_input_sequences", "pct_of_bin", "consensus_length", "gc_content",
               "total_copy_number"]
    rows_out = []
    for i in range(n):
        motif, period, window = args.motifs[i], args.periods[i], args.windows[i]
        for r in read_tsv(args.summaries[i]):
            rows_out.append({
                "label": r["label"], "motif": motif, "target_period": period, "window": window,
                "rank": r["rank"], "n_input_sequences": r["n_input_sequences"],
                "pct_of_bin": r["pct_of_bin"], "consensus_length": r["consensus_length"],
                "gc_content": r["gc_content"], "total_copy_number": r["total_copy_number"],
            })

    rows_out.sort(key=lambda r: int(r["n_input_sequences"]), reverse=True)

    with open(args.out_tsv, "w") as out:
        out.write("\t".join(columns) + "\n")
        for row in rows_out:
            out.write("\t".join(str(row[c]) for c in columns) + "\n")

    sys.stderr.write(
        f"[build_all_clusters_table] {n} bins scanned, {len(rows_out)} total clusters found\n"
    )


if __name__ == "__main__":
    main()
