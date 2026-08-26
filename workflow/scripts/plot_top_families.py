#!/usr/bin/env python
"""
plot_top_families.py

Stage 6: reads the final summary_table.tsv (resolve_redundancy.py's
output) and produces family-ranking bar charts:
  - one plot per individual, ranking that individual's own families by
    copy number, top --top-n. Since an individual can have multiple
    entries (e.g. a hifiasm entry and a verkko entry), the value plotted
    per family is the MAX source_copy_number across that individual's own
    entries for that family — not a sum, since summing would conflate two
    different methods' independent measurements of what should be the
    same underlying biological quantity.
  - one global plot, same max-based aggregation but across ALL entries
    (all individuals combined), top --top-n.
Bars in both are colored by whether the family was found in more than one
individual (n_individuals_present > 1 in summary_table.tsv), so the plot
doubles as a visual cross-check against that column.
"""
import argparse
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COLOR_MULTI_INDIVIDUAL = "#4C72B0"
COLOR_SINGLE_INDIVIDUAL = "#DD8452"


def read_tsv(path):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = [dict(zip(header, line.rstrip("\n").split("\t"))) for line in f if line.strip()]
    return rows


def plot_ranking(families, title, out_png, top_n):
    """families: list of (final_name, value, n_individuals_present), already
    filtered to the families relevant for this plot."""
    ranked = sorted(families, key=lambda f: f[1], reverse=True)[:top_n]
    ranked = ranked[::-1]  # matplotlib barh draws bottom-to-top; reverse so rank 1 is on top

    fig, ax = plt.subplots(figsize=(8, max(3, 0.35 * len(ranked) + 1)))
    if not ranked:
        ax.text(0.5, 0.5, "no families to plot", ha="center", va="center")
    else:
        colors = [COLOR_MULTI_INDIVIDUAL if n > 1 else COLOR_SINGLE_INDIVIDUAL for _, _, n in ranked]
        ax.barh([f[0] for f in ranked], [f[1] for f in ranked], color=colors)
        ax.set_xlabel("copy number (max across entries)")
        legend_handles = [
            plt.Rectangle((0, 0), 1, 1, color=COLOR_MULTI_INDIVIDUAL, label="found in >1 individual"),
            plt.Rectangle((0, 0), 1, 1, color=COLOR_SINGLE_INDIVIDUAL, label="found in 1 individual"),
        ]
        ax.legend(handles=legend_handles, loc="lower right")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--summary-tsv", required=True, help="summary_table.tsv")
    ap.add_argument("--individuals", nargs="+", required=True)
    ap.add_argument("--out-per-individual", nargs="+", required=True,
                     help="Output PNG paths, positionally paired with --individuals")
    ap.add_argument("--out-global", required=True)
    ap.add_argument("--top-n", type=int, default=20)
    args = ap.parse_args()

    if len(args.out_per_individual) != len(args.individuals):
        sys.exit("ERROR: --out-per-individual must have one path per --individuals entry")

    rows = read_tsv(args.summary_tsv)

    n_individuals_present = {}
    per_individual_values = defaultdict(lambda: defaultdict(list))  # individual -> final_name -> [copy_numbers]
    global_values = defaultdict(list)  # final_name -> [copy_numbers]

    for r in rows:
        final_name = r["final_name"]
        copy_number = float(r["source_copy_number"])
        n_individuals_present[final_name] = int(r["n_individuals_present"])
        per_individual_values[r["source_individual"]][final_name].append(copy_number)
        global_values[final_name].append(copy_number)

    for individual, out_png in zip(args.individuals, args.out_per_individual):
        families = [
            (final_name, max(values), n_individuals_present[final_name])
            for final_name, values in per_individual_values.get(individual, {}).items()
        ]
        plot_ranking(families, f"Top families — {individual}", out_png, args.top_n)

    global_families = [
        (final_name, max(values), n_individuals_present[final_name])
        for final_name, values in global_values.items()
    ]
    plot_ranking(global_families, "Top families — all individuals/entries", args.out_global, args.top_n)

    sys.stderr.write(
        f"[plot_top_families] {len(global_values)} families total, plotted top "
        f"{min(args.top_n, len(global_values))} globally and per-individual for "
        f"{len(args.individuals)} individual(s)\n"
    )


if __name__ == "__main__":
    main()
