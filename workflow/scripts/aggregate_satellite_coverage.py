#!/usr/bin/env python
"""
Combines every sample's coverage_summary.tsv into one table, plus a
genome-wide histogram of likely-satellite region sizes split by covered vs.
not -- the "did we lose it somewhere along the way" plot.
"""
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_tsv(path):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in f]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--summaries", nargs="+", required=True)
    ap.add_argument("--regions", nargs="+", required=True)
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--out-histogram", required=True)
    args = ap.parse_args()

    all_summaries = [r for p in args.summaries for r in read_tsv(p)]
    total_bp = sum(int(r["total_likely_bp"]) for r in all_summaries)
    covered_bp = sum(int(r["covered_bp"]) for r in all_summaries)
    n_regions = sum(int(r["n_regions"]) for r in all_summaries)
    n_covered = sum(int(r["n_covered_regions"]) for r in all_summaries)

    cols = ["sample", "total_likely_bp", "covered_bp", "overall_bp_recall", "n_regions", "n_covered_regions", "region_recall"]
    with open(args.out_summary, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in all_summaries:
            f.write("\t".join(r[c] for c in cols) + "\n")
        overall_recall = covered_bp / total_bp if total_bp else 0.0
        region_recall = n_covered / n_regions if n_regions else 0.0
        f.write(f"ALL\t{total_bp}\t{covered_bp}\t{overall_recall:.4f}\t{n_regions}\t{n_covered}\t{region_recall:.4f}\n")

    all_regions = [r for p in args.regions for r in read_tsv(p)]
    lengths_covered = [int(r["length"]) for r in all_regions if r["covered"] == "True"]
    lengths_missing = [int(r["length"]) for r in all_regions if r["covered"] == "False"]

    fig, ax = plt.subplots(figsize=(8, 5))
    all_lengths = lengths_covered + lengths_missing
    bins = np.logspace(np.log10(max(min(all_lengths), 1)), np.log10(max(all_lengths)), 30) if all_lengths else 30
    ax.hist([lengths_covered, lengths_missing], bins=bins, stacked=True,
            label=[f"covered (n={len(lengths_covered)})", f"not covered (n={len(lengths_missing)})"],
            color=["#4C72B0", "#C44E52"])
    ax.set_xscale("log")
    ax.set_xlabel("Likely-satellite region size (bp)")
    ax.set_ylabel("Number of regions")
    ax.set_title("TRF-flagged satellite regions: covered by discovered library vs. not")
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.out_histogram, dpi=150)


if __name__ == "__main__":
    main()
