#!/usr/bin/env python
"""
Parses a RepeatMasker .out file into BED (1-based inclusive -> 0-based
half-open), recoding the "C" strand column to "-". No class/family
filtering -- the point here is "did any library sequence hit this region,"
not class purity (unlike satellite_arrays' own rm_to_bed.py).
"""
import argparse


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rm-out", required=True)
    ap.add_argument("--out-bed", required=True)
    args = ap.parse_args()

    with open(args.rm_out) as f, open(args.out_bed, "w") as out:
        for line in f:
            fields = line.split()
            if len(fields) < 11 or not fields[0].isdigit():
                continue
            query_seq = fields[4]
            query_begin, query_end = int(fields[5]), int(fields[6])
            strand = "-" if fields[8] == "C" else "+"
            name, family, score = fields[9], fields[10], fields[0]
            out.write(f"{query_seq}\t{query_begin - 1}\t{query_end}\t{name}#{family}\t{score}\t{strand}\n")


if __name__ == "__main__":
    main()
