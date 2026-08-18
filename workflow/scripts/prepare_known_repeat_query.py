#!/usr/bin/env python
"""
prepare_known_repeat_query.py

Produce a RepeatMasker-safe query FASTA from repeatmasker_custom_lib.fasta
(this pipeline's final output — every cluster, winners and redundancy-
demoted ones both) for the known-repeat screening step.

Each header in repeatmasker_custom_lib.fasta is ">{label}#{classification}"
(or ">{label}#{classification}/redundant_with_{winner}" for demoted
clusters) — everything from the first "#" onward is our own classification
metadata, not part of the label, and risky to leave in as-is since "#" has
special meaning in RepeatMasker's own library format. This script strips it
down to the bare label.

With --double, each record's sequence is written as seq+seq instead of seq,
so a known-repeat database entry that starts at a different rotation phase
than our consensus can still align end-to-end against a contiguous stretch
of the query, rather than being missed because the match wraps around the
end of a single, un-doubled copy. Same trick used internally elsewhere in
this pipeline (see cross_motif_comparison.py).
"""
import argparse


def read_fasta(path):
    ids, seqs = [], []
    seq = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if ids:
                    seqs.append("".join(seq))
                ids.append(line[1:].strip())
                seq = []
            else:
                seq.append(line.strip())
        if ids:
            seqs.append("".join(seq))
    return list(zip(ids, seqs))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in-fasta", required=True, help="results/repeatmasker_custom_lib.fasta")
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--double", action="store_true",
                     help="Write seq+seq instead of seq for each record")
    args = ap.parse_args()

    records = read_fasta(args.in_fasta)

    with open(args.out_fasta, "w") as out:
        for header, seq in records:
            label = header.split("#", 1)[0].strip()
            out_seq = seq + seq if args.double else seq
            out.write(f">{label}\n{out_seq}\n")


if __name__ == "__main__":
    main()
