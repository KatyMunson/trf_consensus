#!/usr/bin/env python
"""
plot_copy_number_vs_recurrence.py

Runs after Stage 5 (resolve_redundancy.py's summary_table.tsv exists).
Purpose: a family with a middling copy number that nonetheless turns up
independently across several entries/methods is much stronger evidence of
a real satellite than an equally-scored family found only once --
cross-entry recurrence is a form of confirmation a single-entry
copy-number threshold can't see at all. This is what makes a more
permissive min_total_copy_number a defensible choice in combination with
this diagnostic: rather than asking one static per-entry threshold to do
all the noise/signal separation on its own, borderline-copy-number
clusters get a second chance to prove themselves via recurrence, and this
plot checks whether that's actually happening on real data.

One point per family (the representative row, is_representative == True),
total_copy_number (that representative's own source_copy_number, log
x-axis) vs. recurrence -- shown as two panels sharing that x-axis, since
n_entries_found and n_methods_confirming can tell different stories (e.g.
a family found in 4 entries but only 1 method is a weaker claim than one
found in 2 entries spanning 2 different methods).

No new filtering or ranking logic -- purely descriptive, reads
summary_table.tsv as-is.
"""
import argparse
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MARKER_COLOR = "#4C72B0"


def read_tsv(path):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = [dict(zip(header, line.rstrip("\n").split("\t"))) for line in f if line.strip()]
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--summary-tsv", required=True, help="summary_table.tsv")
    ap.add_argument("--out-png", required=True)
    args = ap.parse_args()

    rows = read_tsv(args.summary_tsv)
    representatives = [r for r in rows if r["is_representative"] == "True"]

    if not representatives:
        sys.exit(f"ERROR: no is_representative=True rows found in {args.summary_tsv} -- nothing to plot")

    copy_number = [float(r["source_copy_number"]) for r in representatives]
    n_entries_found = [int(r["n_entries_found"]) for r in representatives]
    n_methods_confirming = [int(r["n_methods_confirming"]) for r in representatives]

    fig, axes = plt.subplots(1, 2, figsize=(11, 5.5))
    fig.suptitle("Post-harmonization: family copy number vs. cross-entry recurrence", y=0.98)

    panels = [
        (axes[0], n_entries_found, "n_entries_found", "copy number vs. entries confirming"),
        (axes[1], n_methods_confirming, "n_methods_confirming", "copy number vs. methods confirming"),
    ]
    for ax, y_values, y_label, title in panels:
        ax.scatter(copy_number, y_values, c=MARKER_COLOR, alpha=0.7, edgecolors="none")
        ax.set_xscale("log")
        ax.set_xlabel("total_copy_number (representative, log)")
        ax.set_ylabel(y_label)
        ax.set_title(title)
        max_y = max(y_values)
        ax.set_yticks(range(1, max_y + 1))

    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(args.out_png, dpi=150)

    sys.stderr.write(
        f"[plot_copy_number_vs_recurrence] {len(representatives)} families plotted, "
        f"wrote {args.out_png}\n"
    )


if __name__ == "__main__":
    main()
