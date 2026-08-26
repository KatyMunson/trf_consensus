#!/usr/bin/env python
"""
BP-level overlap between likely-satellite regions and this sample's -lib
hits. A region counts as "covered" if its overlap fraction >=
--min-coverage-fraction (this classification is only used for the
per-region covered/not-covered breakdown and the missing-regions list --
the overall bp recall in coverage_summary.tsv sums actual overlapping bp
across all regions regardless of the per-region threshold).
"""
import argparse
from collections import defaultdict


def read_bed(path):
    by_seq = defaultdict(list)
    with open(path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            by_seq[fields[0]].append((int(fields[1]), int(fields[2])))
    for seq in by_seq:
        by_seq[seq].sort()
    return by_seq


def merged_overlap_bp(region_start, region_end, hit_intervals):
    covered = 0
    cur_end = region_start
    for h_start, h_end in hit_intervals:
        if h_end <= region_start or h_start >= region_end:
            continue
        ov_start = max(h_start, region_start, cur_end)
        ov_end = min(h_end, region_end)
        if ov_end > ov_start:
            covered += ov_end - ov_start
            cur_end = max(cur_end, ov_end)
    return covered


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--likely-bed", required=True)
    ap.add_argument("--hits-bed", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--min-coverage-fraction", type=float, required=True)
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--out-missing", required=True)
    ap.add_argument("--out-regions", required=True)
    args = ap.parse_args()

    likely = read_bed(args.likely_bed)
    hits = read_bed(args.hits_bed)

    total_bp = covered_bp = n_regions = n_covered_regions = 0

    with open(args.out_regions, "w") as regions_out, open(args.out_missing, "w") as missing_out:
        regions_out.write("sample\tseqname\tstart\tend\tlength\tcovered_bp\tcoverage_fraction\tcovered\n")
        missing_out.write("sample\tseqname\tstart\tend\tlength\tcovered_bp\tcoverage_fraction\n")

        for seqname, ivs in sorted(likely.items()):
            hit_ivs = hits.get(seqname, [])
            for start, end in ivs:
                length = end - start
                cov_bp = merged_overlap_bp(start, end, hit_ivs)
                frac = cov_bp / length if length else 0.0
                is_covered = frac >= args.min_coverage_fraction

                total_bp += length
                covered_bp += cov_bp
                n_regions += 1
                n_covered_regions += int(is_covered)

                regions_out.write(f"{args.sample}\t{seqname}\t{start}\t{end}\t{length}\t{cov_bp}\t{frac:.4f}\t{is_covered}\n")
                if not is_covered:
                    missing_out.write(f"{args.sample}\t{seqname}\t{start}\t{end}\t{length}\t{cov_bp}\t{frac:.4f}\n")

    overall_recall = covered_bp / total_bp if total_bp else 0.0
    region_recall = n_covered_regions / n_regions if n_regions else 0.0
    with open(args.out_summary, "w") as f:
        f.write("sample\ttotal_likely_bp\tcovered_bp\toverall_bp_recall\tn_regions\tn_covered_regions\tregion_recall\n")
        f.write(f"{args.sample}\t{total_bp}\t{covered_bp}\t{overall_recall:.4f}\t{n_regions}\t{n_covered_regions}\t{region_recall:.4f}\n")


if __name__ == "__main__":
    main()
