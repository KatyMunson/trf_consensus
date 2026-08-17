#!/usr/bin/env python
"""
extract_by_period.py

Pull entries from a native TRF .dat file whose period size (field 3) falls
within [period - window, period + window], and write the per-entry consensus
motif (field 14) to a FASTA file. A companion TSV records provenance/QC
metadata (source sequence, coordinates, copy number, consensus size, percent
match) so entries can be traced back to the .dat file later.

TRF .dat data lines (whitespace-delimited):
  1 start, 2 end, 3 period_size, 4 copy_number, 5 consensus_size,
  6 percent_matches, 7 percent_indels, 8 alignment_score,
  9-12 %A %C %G %T, 13 entropy, 14 consensus_sequence,
  15 repeat_sequence, [16 left_flank, 17 right_flank]
Header lines are used only to track which contig subsequent entries belong
to, and are recognized in two forms:
  - Standard TRF: "Sequence: <name>" (+ a separate "Parameters:" line)
  - TRF -ngs mode: "@<name>" as its own line, no "Parameters:" line
Neither form is whitespace-delimited into >=15 numeric-leading fields, so
both are naturally skipped as data rows regardless. If the .dat file has
neither header form at all (e.g. headers were stripped when it was
generated/concatenated), pass --default-seqname so entries are labeled with
the real contig/scaffold name instead of "unknown".
"""
import argparse
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dat", required=True, help="Input TRF .dat file")
    p.add_argument("--period", type=int, required=True, help="Target period size")
    p.add_argument("--window", type=int, required=True, help="+/- tolerance on period")
    p.add_argument("--out-fasta", required=True, help="Output FASTA of consensus motifs")
    p.add_argument("--out-tsv", required=True, help="Output TSV of provenance/QC metadata")
    p.add_argument("--default-seqname", default=None,
                    help="Source sequence name to use for entries seen before any "
                         "header line, or for the whole file if it has none at all. "
                         "Leave unset to fall back to 'unknown'.")
    return p.parse_args()


def main():
    args = parse_args()
    lo, hi = args.period - args.window, args.period + args.window

    current_seq = args.default_seqname if args.default_seqname else "unknown"
    saw_header = False
    n_written = 0
    n_lines = 0

    with open(args.dat) as fin, open(args.out_fasta, "w") as fa, open(args.out_tsv, "w") as tsv:
        tsv.write("id\tsource_seq\tstart\tend\tperiod\tcopy_number\tconsensus_size\tpercent_match\n")

        for lineno, line in enumerate(fin, 1):
            line = line.rstrip("\n")
            if not line.strip():
                continue

            if line.startswith("Sequence:"):
                # Standard TRF: "Sequence: h2tg000001l"
                parts = line.split(None, 1)
                current_seq = parts[1].strip() if len(parts) > 1 else "unknown"
                saw_header = True
                continue

            if line.startswith("@"):
                # TRF -ngs mode: "@h2tg000001l"
                tokens = line[1:].split()
                current_seq = tokens[0] if tokens else "unknown"
                saw_header = True
                continue

            if line.startswith("Parameters:"):
                continue

            fields = line.split()
            # Data lines start with two integers (start, end); guard against
            # any other stray header/annotation lines TRF may emit.
            if len(fields) < 15 or not (fields[0].isdigit() and fields[1].isdigit()):
                continue

            n_lines += 1
            try:
                start, end = fields[0], fields[1]
                period = int(round(float(fields[2])))
                copy_number = fields[3]
                consensus_size = fields[4]
                percent_match = fields[5]
                consensus_seq = fields[13]
            except (IndexError, ValueError):
                sys.stderr.write(f"WARNING: could not parse .dat line {lineno}, skipping\n")
                continue

            if lo <= period <= hi:
                entry_id = f"{current_seq}:{start}-{end}"
                fa.write(f">{entry_id}\n{consensus_seq}\n")
                tsv.write(
                    f"{entry_id}\t{current_seq}\t{start}\t{end}\t{period}\t"
                    f"{copy_number}\t{consensus_size}\t{percent_match}\n"
                )
                n_written += 1

    sys.stderr.write(
        f"[extract_by_period] scanned {n_lines} data lines, "
        f"wrote {n_written} entries with period in [{lo},{hi}]\n"
    )
    if n_written == 0:
        sys.stderr.write("WARNING: no entries matched — check --period/--window and .dat field layout\n")
    if not saw_header and current_seq == "unknown":
        sys.stderr.write(
            "WARNING: no 'Sequence:' header lines found in .dat and no --default-seqname "
            "given — every entry is labeled source_seq='unknown'. Pass --default-seqname "
            "if you know the contig/scaffold this .dat came from.\n"
        )


if __name__ == "__main__":
    main()
