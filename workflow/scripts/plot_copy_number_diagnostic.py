#!/usr/bin/env python
"""
plot_copy_number_diagnostic.py

Stage 2b: runs once every manifest entry has completed Stage 2
(build_all_clusters_table.py), before the Stage 3 filter is applied. Pools
every entry's total_copy_number values -- pre-filter, including clusters
that are about to be dropped, that's the whole point -- into one
log-x-scaled histogram, with a vertical reference line at the configured
min_total_copy_number. Purpose: eyeball whether the threshold actually
sits in a real gap between background noise and true arrays before
committing to a rerun.
"""
import argparse
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_tsv(path):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = [dict(zip(header, line.rstrip("\n").split("\t"))) for line in f if line.strip()]
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--clusters-tsv", nargs="+", required=True,
                     help="Every entry's all_clusters.tsv (pre-filter)")
    ap.add_argument("--min-total-copy-number", type=float, required=True)
    ap.add_argument("--out-png", required=True)
    args = ap.parse_args()

    values = []
    for path in args.clusters_tsv:
        for row in read_tsv(path):
            values.append(float(row["total_copy_number"]))

    if not values:
        sys.exit("ERROR: no clusters found across --clusters-tsv inputs -- nothing to plot")

    values = np.array(values)
    if (values <= 0).any():
        sys.exit("ERROR: total_copy_number must be strictly positive for a log-scaled histogram")

    fig, ax = plt.subplots(figsize=(8, 5))
    bins = np.logspace(np.log10(values.min()), np.log10(values.max()), 40)
    ax.hist(values, bins=bins, color="#4C72B0", edgecolor="white", linewidth=0.5)
    ax.set_xscale("log")
    ax.axvline(args.min_total_copy_number, color="#C44E52", linestyle="--", linewidth=1.5,
               label=f"min_total_copy_number = {args.min_total_copy_number:g}")
    ax.set_xlabel("total_copy_number (log scale)")
    ax.set_ylabel("number of clusters")
    ax.set_title("Pre-filter cluster copy-number distribution, all entries pooled")
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.out_png, dpi=150)

    n_below = int((values < args.min_total_copy_number).sum())
    sys.stderr.write(
        f"[plot_copy_number_diagnostic] {len(values)} clusters pooled across "
        f"{len(args.clusters_tsv)} entries, {n_below} below threshold "
        f"{args.min_total_copy_number:g}, wrote {args.out_png}\n"
    )


if __name__ == "__main__":
    main()
