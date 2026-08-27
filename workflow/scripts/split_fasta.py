#!/usr/bin/env python3
"""
split_fasta.py -- snake/boustrophedon split a FASTA file's records across N
output files, by contig (record) count, not by base count.

Adapted from compare_meadowlark's scripts/split_fasta.py (itself adapted
from vendor/rhodonite's split_fasta.py), upgraded from plain round-robin
(record i -> output i % N, always wrapping N-1 -> 0) to a snake/boustrophedon
assignment: the target index bounces back and forth (0,1,...,N-1,N-1,...,1,0,
0,1,...) instead of always wrapping. Still a single sequential line-scan, no
sequence-length lookahead, no pysam -- this only actually improves load
balance over plain round-robin when contigs happen to be roughly
size-ordered in the source fasta already (common for hifiasm/verkko output
-- largest contigs first -- but not guaranteed); if the input order is
arbitrary, it's no worse than plain round-robin. Not a substitute for real
bp-balanced bin-packing, which remains out of scope here.
"""
import argparse
import sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--infile", required=True, help="input FASTA file")
    parser.add_argument(
        "--outputs", nargs="+", required=True, help="output FASTA chunk files"
    )
    args = parser.parse_args()

    n = len(args.outputs)
    outs = [open(f, "w") for f in args.outputs]
    current = None
    next_idx = 0
    direction = 1
    with open(args.infile) as fasta:
        for line in fasta:
            if line.startswith(">"):
                current = outs[next_idx]
                next_idx += direction
                if next_idx == n:
                    next_idx = n - 1
                    direction = -1
                elif next_idx < 0:
                    next_idx = 0
                    direction = 1
            if current is None:
                sys.exit(
                    f"ERROR: {args.infile} has sequence data before its "
                    "first '>' header"
                )
            current.write(line)

    for out in outs:
        out.close()
