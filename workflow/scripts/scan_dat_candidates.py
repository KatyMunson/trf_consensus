#!/usr/bin/env python
"""
scan_dat_candidates.py

Scan one or more raw TRF .dat files and nominate candidate period (monomer)
sizes worth building a consensus-of-consensus for.

A candidate period passes if EITHER:
  A. Distributed evidence — at least --min-blocks total loci fall within
     period +/- window, AND at least one of them has copy_number >=
     --min-copy-number
  B. Single massive block — ANY one locus within the window has copy_number
     >= --min-single-block-copy-number, regardless of how many total loci
     are in the window

(B) exists because a genuine large satellite/centromeric array is often
captured by TRF as ONE contiguous locus (e.g. a single 221kb, 534-copy
block) rather than fragmenting into several — requiring --min-blocks would
wrongly reject exactly the kind of candidate this is meant to find.
--min-single-block-copy-number defaults higher than --min-copy-number so a
lone modest-copy-number block doesn't pass on its own.

--min-period-length filters out micro/mini-satellites, keeping only larger
tandem repeat classes (default 150bp; human alpha-satellite is 171bp, but
other characterized satellite families run shorter — consider lowering for
a permissive discovery pass).

Every eligible period value is scored on its own +/-window neighborhood
first (independently — not merged with its neighbors). Candidates are then
selected by non-maximum suppression: take the best-scoring remaining period,
report it, and suppress every other period within +/-nms-radius of it
(rather than chaining merges transitively across a whole dense run of
adjacent period values, which is what a naive "merge if gap<=window" scheme
does — real .dat files often have a near-continuous run of distinct period
values across hundreds of bp, and a chained merge collapses all of them,
including any true candidates in the middle, into one blob keyed on
whichever value happens to have the most raw loci). This keeps candidates
that are >nms-radius apart independent even if the space between them is
densely populated with unrelated periods.

Reported per-candidate stats include total_array_bp (sum of end-start across
all loci in the window — the single most direct "how much genome does this
occupy" number for prioritizing true large-satellite candidates) and
mean_percent_match / mean_entropy (TRF's own alignment-quality and Shannon-
entropy fields — genuine satellite monomers tend to sit well above 0 entropy
on TRF's 0-2 bit scale; homopolymer/microsatellite noise sits much lower).

Three outputs:
  --out-tsv        final candidates after non-maximum suppression
  --out-raw-tsv    (optional) EVERY eligible period's own individual window
                   stats, unsuppressed — use this to look up a specific
                   period you expected to see and check why it either did
                   or didn't survive suppression
  --out-manifest   (optional) passing candidates in periods.tsv format
"""
import argparse
import sys
from collections import defaultdict


def parse_dat_entries(paths):
    """Yield dicts for every data line across all given .dat files. Header
    lines are recognized in both 'Sequence: <name>' and '@<name>' (-ngs
    mode) forms, purely for source_seq bookkeeping."""
    for path in paths:
        current_seq = "unknown"
        with open(path) as fin:
            for line in fin:
                line = line.rstrip("\n")
                if not line.strip():
                    continue
                if line.startswith("Sequence:"):
                    parts = line.split(None, 1)
                    current_seq = parts[1].strip() if len(parts) > 1 else "unknown"
                    continue
                if line.startswith("@"):
                    tokens = line[1:].split()
                    current_seq = tokens[0] if tokens else "unknown"
                    continue
                if line.startswith("Parameters:"):
                    continue

                fields = line.split()
                if len(fields) < 15 or not (fields[0].isdigit() and fields[1].isdigit()):
                    continue
                try:
                    yield {
                        "source_seq": current_seq,
                        "start": int(fields[0]),
                        "end": int(fields[1]),
                        "period": int(round(float(fields[2]))),
                        "copy_number": float(fields[3]),
                        "percent_match": float(fields[5]),
                        "entropy": float(fields[12]),
                    }
                except (IndexError, ValueError):
                    continue


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dat", nargs="+", required=True, help="One or more TRF .dat files")
    ap.add_argument("--min-copy-number", type=float, default=100.0,
                     help="Distributed-evidence rule: a candidate needs >=1 window locus with "
                          "copy_number at or above this")
    ap.add_argument("--min-blocks", type=int, default=5,
                     help="Distributed-evidence rule: a candidate needs >= this many total loci "
                          "within period +/- window")
    ap.add_argument("--min-single-block-copy-number", type=float, default=300.0,
                     help="Single-massive-block rule: a candidate passes if ANY one locus in the "
                          "window has copy_number at or above this, regardless of block count. "
                          "Set higher than --min-copy-number so a lone modest block doesn't pass "
                          "on its own; set to a very large number (e.g. 1e9) to disable this rule "
                          "and require distributed evidence only.")
    ap.add_argument("--min-period-length", type=int, default=150,
                     help="Candidate period itself must be >= this many bp")
    ap.add_argument("--window", type=int, default=5,
                     help="+/- tolerance for counting blocks per candidate period; "
                          "matches the main pipeline's window")
    ap.add_argument("--nms-radius", type=int, default=None,
                     help="After picking a candidate, suppress other period values within "
                          "+/- this many bp of it. Defaults to --window if not set. Keep this "
                          "well below the typical spacing between your real repeat families, "
                          "or true candidates closer together than this will suppress each other.")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-raw-tsv", default=None,
                     help="Optional: every eligible period's own stats, before suppression")
    ap.add_argument("--out-manifest", default=None,
                     help="Optional: also write candidates in periods.tsv format (name, period, window)")
    ap.add_argument("--manifest-name-prefix", default="cand",
                     help="Prefix for auto-generated names in --out-manifest, e.g. cand412")
    args = ap.parse_args()

    nms_radius = args.nms_radius if args.nms_radius is not None else args.window

    entries = list(parse_dat_entries(args.dat))
    sys.stderr.write(f"[scan_dat_candidates] parsed {len(entries)} data lines from {len(args.dat)} file(s)\n")

    by_period = defaultdict(list)
    for e in entries:
        by_period[e["period"]].append(e)

    eligible_periods = sorted(p for p in by_period if p >= args.min_period_length)

    raw_columns = ["period", "window", "n_blocks", "max_copy_number",
                   "n_loci_ge_min_copy_number", "n_distinct_source_seqs",
                   "mean_array_length_bp", "total_array_bp", "mean_percent_match",
                   "mean_entropy", "passes_distributed", "passes_single_block", "passes_filters"]
    final_columns = ["candidate_period", "window", "n_blocks", "max_copy_number",
                      "n_loci_ge_min_copy_number", "n_distinct_source_seqs",
                      "mean_array_length_bp", "total_array_bp", "mean_percent_match",
                      "mean_entropy", "suppressed_period_range",
                      "n_periods_suppressed", "passes_distributed", "passes_single_block",
                      "passes_filters"]

    if not eligible_periods:
        sys.stderr.write("WARNING: no periods meet --min-period-length; writing empty outputs\n")
        with open(args.out_tsv, "w") as out:
            out.write("\t".join(final_columns) + "\n")
        if args.out_raw_tsv:
            with open(args.out_raw_tsv, "w") as out:
                out.write("\t".join(raw_columns) + "\n")
        if args.out_manifest:
            open(args.out_manifest, "w").close()
        return

    # --- score every eligible period independently on its own +/-window ---
    period_stats = {}
    for p in eligible_periods:
        lo, hi = p - args.window, p + args.window
        window_entries = [e for q in range(lo, hi + 1) for e in by_period.get(q, [])]
        n_blocks = len(window_entries)
        if n_blocks == 0:
            continue
        max_copy = max(e["copy_number"] for e in window_entries)
        n_ge_min_copy = sum(1 for e in window_entries if e["copy_number"] >= args.min_copy_number)
        n_source_seqs = len({e["source_seq"] for e in window_entries})
        total_bp = sum(e["end"] - e["start"] for e in window_entries)
        mean_len = total_bp / n_blocks
        mean_pct_match = sum(e["percent_match"] for e in window_entries) / n_blocks
        mean_entropy = sum(e["entropy"] for e in window_entries) / n_blocks

        passes_distributed = (max_copy >= args.min_copy_number) and (n_blocks >= args.min_blocks)
        passes_single_block = max_copy >= args.min_single_block_copy_number
        period_stats[p] = {
            "n_blocks": n_blocks,
            "max_copy_number": max_copy,
            "n_loci_ge_min_copy_number": n_ge_min_copy,
            "n_distinct_source_seqs": n_source_seqs,
            "mean_array_length_bp": round(mean_len, 1),
            "total_array_bp": total_bp,
            "mean_percent_match": round(mean_pct_match, 1),
            "mean_entropy": round(mean_entropy, 2),
            "passes_distributed": passes_distributed,
            "passes_single_block": passes_single_block,
            "passes_filters": passes_distributed or passes_single_block,
        }

    if args.out_raw_tsv:
        with open(args.out_raw_tsv, "w") as out:
            out.write("\t".join(raw_columns) + "\n")
            for p in sorted(period_stats):
                s = period_stats[p]
                out.write(f"{p}\t{args.window}\t{s['n_blocks']}\t{s['max_copy_number']}\t"
                          f"{s['n_loci_ge_min_copy_number']}\t{s['n_distinct_source_seqs']}\t"
                          f"{s['mean_array_length_bp']}\t{s['total_array_bp']}\t"
                          f"{s['mean_percent_match']}\t{s['mean_entropy']}\t"
                          f"{s['passes_distributed']}\t{s['passes_single_block']}\t{s['passes_filters']}\n")

    # --- non-maximum suppression: best-quality periods first, suppress only their own vicinity ---
    remaining = set(period_stats)
    ranked = sorted(
        remaining,
        key=lambda p: (period_stats[p]["passes_single_block"],
                        period_stats[p]["n_loci_ge_min_copy_number"],
                        period_stats[p]["total_array_bp"],
                        period_stats[p]["max_copy_number"],
                        -p),
        reverse=True,
    )

    selected = []
    suppressed_by = defaultdict(list)
    for p in ranked:
        if p not in remaining:
            continue
        selected.append(p)
        remaining.discard(p)
        for offset in range(-nms_radius, nms_radius + 1):
            q = p + offset
            if offset != 0 and q in remaining:
                remaining.discard(q)
                suppressed_by[p].append(q)

    selected.sort()
    candidates = []
    for p in selected:
        s = period_stats[p]
        group = sorted([p] + suppressed_by[p])
        candidates.append({
            "candidate_period": p,
            "window": args.window,
            "n_blocks": s["n_blocks"],
            "max_copy_number": s["max_copy_number"],
            "n_loci_ge_min_copy_number": s["n_loci_ge_min_copy_number"],
            "n_distinct_source_seqs": s["n_distinct_source_seqs"],
            "mean_array_length_bp": s["mean_array_length_bp"],
            "total_array_bp": s["total_array_bp"],
            "mean_percent_match": s["mean_percent_match"],
            "mean_entropy": s["mean_entropy"],
            "suppressed_period_range": f"{group[0]}-{group[-1]}",
            "n_periods_suppressed": len(suppressed_by[p]),
            "passes_distributed": s["passes_distributed"],
            "passes_single_block": s["passes_single_block"],
            "passes_filters": s["passes_filters"],
        })

    passing = [c for c in candidates if c["passes_filters"]]

    with open(args.out_tsv, "w") as out:
        out.write("\t".join(final_columns) + "\n")
        for c in candidates:
            out.write("\t".join(str(c[col]) for col in final_columns) + "\n")

    if args.out_manifest:
        with open(args.out_manifest, "w") as out:
            out.write("# name\tperiod\twindow\n")
            for c in passing:
                name = f"{args.manifest_name_prefix}{c['candidate_period']}"
                out.write(f"{name}\t{c['candidate_period']}\t{c['window']}\n")

    sys.stderr.write(
        f"[scan_dat_candidates] {len(eligible_periods)} distinct periods >= "
        f"{args.min_period_length}bp scored individually, {len(candidates)} candidates "
        f"survived non-maximum suppression (radius={nms_radius}), "
        f"{len(passing)} pass (distributed: min_copy_number>={args.min_copy_number} and "
        f"min_blocks>={args.min_blocks}) OR (single block: copy_number>={args.min_single_block_copy_number})\n"
    )


if __name__ == "__main__":
    main()
