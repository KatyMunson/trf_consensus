#!/usr/bin/env python
"""
plot_copy_number_qc_diagnostic.py

Runs alongside Stage 2b's plot_copy_number_diagnostic.py (same input
timing: after every entry's all_clusters.tsv exists, before the Stage 3
filter is applied). The plain total_copy_number histogram alone can't
resolve an ambiguous middle -- a noisy plateau with no clean valley, say --
so this cross-references the same pre-filter total_copy_number against
three other axes: target_period (already available per cluster) and two
of TRF's own per-locus quality fields, aggregated per cluster by
rank_family_clusters.py (mean_percent_match, mean_entropy, both
copy-number-weighted -- see that script's docstring).

Three scatter panels, one row, sharing a log-scaled total_copy_number
x-axis. Points are colored by whether they currently pass
min_total_copy_number (blue) or not (red) -- this isn't a new
classification, just today's threshold decision made visible against the
other axes, so disagreements are easy to spot:
  - red points sitting high on percent-match/entropy: candidates the
    current threshold may be wrongly excluding (e.g. an array fragmented
    across assembly gaps -- low copy number for a reason other than noise)
  - blue points sitting low on percent-match/entropy: candidates the
    threshold is letting through that look questionable on every other axis
"""
import argparse
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COLOR_PASS = "#4C72B0"
COLOR_FAIL = "#C44E52"


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

    total_copy_number, target_period, mean_percent_match, mean_entropy = [], [], [], []
    for path in args.clusters_tsv:
        for row in read_tsv(path):
            total_copy_number.append(float(row["total_copy_number"]))
            target_period.append(float(row["target_period"]))
            mean_percent_match.append(float(row["mean_percent_match"]))
            mean_entropy.append(float(row["mean_entropy"]))

    if not total_copy_number:
        sys.exit("ERROR: no clusters found across --clusters-tsv inputs -- nothing to plot")

    passes = [cn >= args.min_total_copy_number for cn in total_copy_number]
    colors = [COLOR_PASS if p else COLOR_FAIL for p in passes]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5.5))
    fig.suptitle("Pre-filter QC cross-reference — copy number alone vs. TRF's own quality fields",
                 y=0.98)

    panels = [
        (axes[0], target_period, "target_period", "copy number vs. period length"),
        (axes[1], mean_percent_match, "mean_percent_match", "copy number vs. percent match"),
        (axes[2], mean_entropy, "mean_entropy", "copy number vs. entropy"),
    ]
    for ax, y_values, y_label, title in panels:
        ax.scatter(total_copy_number, y_values, c=colors, alpha=0.7, edgecolors="none")
        ax.set_xscale("log")
        ax.axvline(args.min_total_copy_number, color="black", linestyle="--", linewidth=1)
        ax.set_xlabel("total_copy_number (log)")
        ax.set_ylabel(y_label)
        ax.set_title(title)

    legend_handles = [
        plt.Rectangle((0, 0), 1, 1, color=COLOR_PASS,
                       label=f"total_copy_number ≥ {args.min_total_copy_number:g} (passes today)"),
        plt.Rectangle((0, 0), 1, 1, color=COLOR_FAIL,
                       label=f"total_copy_number < {args.min_total_copy_number:g} (filtered today)"),
    ]
    fig.legend(handles=legend_handles, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 0.90))
    fig.tight_layout(rect=(0, 0, 1, 0.82))
    fig.savefig(args.out_png, dpi=150)

    n_pass = sum(passes)
    sys.stderr.write(
        f"[plot_copy_number_qc_diagnostic] {len(total_copy_number)} clusters pooled across "
        f"{len(args.clusters_tsv)} entries, {n_pass} pass threshold {args.min_total_copy_number:g}, "
        f"{len(total_copy_number) - n_pass} filtered, wrote {args.out_png}\n"
    )


if __name__ == "__main__":
    main()
