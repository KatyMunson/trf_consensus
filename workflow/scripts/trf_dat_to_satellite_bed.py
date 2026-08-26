#!/usr/bin/env python
"""
Flags individual TRF loci that look satellite-like using candidate_scan's
locus-level thresholds (min_copy_number / min_single_block_copy_number /
min_period_length) -- NOT its min_blocks requirement, which nominates
period bins for consensus-building, a different question from whether one
genomic interval looks like satellite DNA. Passing loci within
merge_distance bp on the same sequence are merged into regions.

.dat format (TRF run with -d): "Sequence: <name>" (or "@<name>") header
lines precede each sequence's locus lines. Locus columns:
  start end period_size copy_number consensus_size percent_matches
  percent_indels score pct_A pct_C pct_G pct_T entropy consensus repeat
Only the first 4 numeric columns are used.

Seqname parsing takes only the first whitespace-delimited token after
"Sequence:"/"@", NOT scan_dat_candidates.py's whole-line behavior:
scan_dat_candidates.py only ever uses source_seq for a distinct-count
bookkeeping column, so it doesn't matter there if a header line carries
extra text after the contig name. Here the seqname becomes a BED chrom
field that must match RepeatMasker's query_seq column (which is always
first-token-only), so first-token parsing is the correct behavior for
this script's purpose even though it diverges from scan_dat_candidates.py.
Numeric column indexing (start/end/period/copy_number) is otherwise the
same convention scan_dat_candidates.py uses.
"""
import argparse


def parse_dat(path, default_seqname):
    loci = []
    seqname = default_seqname
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if line.startswith("Sequence:"):
                seqname = line.split("Sequence:", 1)[1].strip().split()[0]
                continue
            if line.startswith("@"):
                seqname = line[1:].strip().split()[0]
                continue
            fields = line.split()
            if len(fields) < 4 or not fields[0].isdigit():
                continue
            start, end = int(fields[0]), int(fields[1])
            period = int(fields[2])
            copy_number = float(fields[3])
            if not seqname:
                raise ValueError(
                    "Locus line with no seqname header and no --default-seqname "
                    "given -- set trf_dat_default_seqname in config.yaml."
                )
            loci.append((seqname, start, end, period, copy_number))
    return loci


def merge_intervals(intervals, merge_distance):
    by_seq = {}
    for seqname, start, end in intervals:
        by_seq.setdefault(seqname, []).append((start, end))
    merged = []
    for seqname, ivs in by_seq.items():
        ivs.sort()
        cur_start, cur_end = ivs[0]
        for start, end in ivs[1:]:
            if start - cur_end <= merge_distance:
                cur_end = max(cur_end, end)
            else:
                merged.append((seqname, cur_start, cur_end))
                cur_start, cur_end = start, end
        merged.append((seqname, cur_start, cur_end))
    return merged


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dat", required=True)
    ap.add_argument("--out-bed", required=True)
    ap.add_argument("--min-copy-number", type=float, required=True)
    ap.add_argument("--min-single-block-copy-number", type=float, required=True)
    ap.add_argument("--min-period-length", type=int, required=True)
    ap.add_argument("--merge-distance", type=int, required=True)
    ap.add_argument("--default-seqname", default="")
    args = ap.parse_args()

    loci = parse_dat(args.dat, args.default_seqname)
    passing = [
        (seqname, start, end)
        for seqname, start, end, period, copy_number in loci
        if period >= args.min_period_length
        and (copy_number >= args.min_copy_number
             or copy_number >= args.min_single_block_copy_number)
    ]
    merged = sorted(merge_intervals(passing, args.merge_distance))

    with open(args.out_bed, "w") as out:
        for seqname, start, end in merged:
            # TRF coords 1-based inclusive -> BED 0-based half-open
            out.write(f"{seqname}\t{start - 1}\t{end}\n")


if __name__ == "__main__":
    main()
