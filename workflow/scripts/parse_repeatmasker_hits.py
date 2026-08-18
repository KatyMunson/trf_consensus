#!/usr/bin/env python
"""
parse_repeatmasker_hits.py

Turn RepeatMasker's .out file (from the known-repeat screening step, run
against our own discovered consensus motifs) into a clean per-label TSV.

Every label in --query-fasta appears as exactly one output row, whether or
not RepeatMasker found anything for it — RepeatMasker's .out file only ever
lists sequences with at least one hit, so the full label universe has to
come from the query FASTA itself, not from the .out file. A missing .out
file (RepeatMasker can decline to create one at all when nothing anywhere
in the query matches anything in the library) or a header-only one is
treated exactly like zero hits found, not an error: for satellite DNA in a
non-model species this is a completely plausible, even expected, real
result.

Parsing approach: RepeatMasker's .out format has a few header/blank lines
before data rows begin, and the exact header formatting has drifted
slightly across versions. Rather than hardcoding a line count, a line is
treated as a data row if, after stripping and splitting on whitespace, its
first token parses as an integer (the SW score column) — the same heuristic
already used in extract_by_period.py for skipping TRF .dat header lines.

Column layout (whitespace-delimited; RepeatMasker's repeat-position columns
near the end swap order depending on strand and aren't needed here):
  [0] SW score            [6] query end
  [1] % divergence        [7] query (left), parenthesized
  [2] % deleted            [8] strand (+ or C)
  [3] % inserted           [9] matching repeat name
  [4] query sequence name  [10] repeat class/family
  [5] query begin          ... (repeat-position fields, not needed)
Strand is normalized from RepeatMasker's "C" to "-" for consistency with
this pipeline's own strand convention elsewhere ("+"/"-").

Aggregation: queries are doubled (seq+seq) when --doubled is set (matching
known_repeat_screen.double_sequences), so up to 2 hit rows can come back
per true label for what's really the same underlying match — one from each
copy in the doubled query. These collapse to one output row per label by
keeping the hit with the highest SW score.

CAVEAT, not something this script tries to correct for: if a real match
happens to span the artificial junction introduced by doubling, RepeatMasker
can report an inflated or split score right at that boundary. This is a
known, accepted approximation of the doubling trick. If a top hit looks
borderline or surprising, it's worth a manual look at the un-doubled
first-half coordinates (query_start/query_end here are reported exactly as
RepeatMasker output them, i.e. in the doubled query's coordinate space when
--doubled was used).
"""
import argparse
import os
import sys

STRAND_MAP = {"+": "+", "C": "-"}


def read_fasta_labels(path):
    """Just the header labels, in file order — this is the complete label
    universe the output table must cover, whether or not each one got a
    RepeatMasker hit."""
    labels = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                labels.append(line[1:].strip().split()[0])
    return labels


def parse_rm_out(path):
    """Yield one dict per data row. Missing file or a file with no data
    rows (header-only, or truly empty) both yield nothing — treated as
    zero hits, not an error."""
    if not path or not os.path.isfile(path):
        sys.stderr.write(
            f"[parse_repeatmasker_hits] no .out file at {path!r} — treating as zero hits "
            f"(RepeatMasker can decline to create one when nothing matches)\n"
        )
        return
    with open(path) as f:
        for line in f:
            fields = line.strip().split()
            if len(fields) < 11 or not fields[0].isdigit():
                continue
            try:
                yield {
                    "sw_score": int(fields[0]),
                    "pct_divergence": fields[1],
                    "pct_deleted": fields[2],
                    "pct_inserted": fields[3],
                    "query_label": fields[4],
                    "query_start": fields[5],
                    "query_end": fields[6],
                    "strand": STRAND_MAP.get(fields[8], fields[8]),
                    "repeat_name": fields[9],
                    "repeat_class_family": fields[10],
                }
            except (IndexError, ValueError):
                continue


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rm-out", required=True, help="RepeatMasker .out file (may not exist)")
    ap.add_argument("--query-fasta", required=True,
                     help="The query FASTA that was screened (prepare_known_repeat_query.py's "
                          "output) — source of the full label universe, since RepeatMasker's "
                          ".out only lists sequences with a hit")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--doubled", action="store_true",
                     help="Matches known_repeat_screen.double_sequences — documentation only, "
                          "see the module docstring's coordinate-space caveat")
    args = ap.parse_args()

    labels = read_fasta_labels(args.query_fasta)

    best_hit = {}
    n_rows = 0
    for row in parse_rm_out(args.rm_out):
        n_rows += 1
        label = row["query_label"]
        if label not in best_hit or row["sw_score"] > best_hit[label]["sw_score"]:
            best_hit[label] = row

    columns = ["label", "has_known_hit", "repeat_name", "repeat_class_family", "sw_score",
               "pct_divergence", "pct_deleted", "pct_inserted", "query_start", "query_end", "strand"]

    with open(args.out_tsv, "w") as out:
        out.write("\t".join(columns) + "\n")
        for label in labels:
            hit = best_hit.get(label)
            if hit is None:
                row = {c: "NA" for c in columns}
                row["label"] = label
                row["has_known_hit"] = False
            else:
                row = {
                    "label": label, "has_known_hit": True,
                    "repeat_name": hit["repeat_name"],
                    "repeat_class_family": hit["repeat_class_family"],
                    "sw_score": hit["sw_score"],
                    "pct_divergence": hit["pct_divergence"],
                    "pct_deleted": hit["pct_deleted"],
                    "pct_inserted": hit["pct_inserted"],
                    "query_start": hit["query_start"],
                    "query_end": hit["query_end"],
                    "strand": hit["strand"],
                }
            out.write("\t".join(str(row[c]) for c in columns) + "\n")

    n_hit_labels = len(best_hit)
    sys.stderr.write(
        f"[parse_repeatmasker_hits] {len(labels)} query labels, {n_rows} raw RepeatMasker hit "
        f"row(s) parsed, {n_hit_labels} label(s) with >=1 hit collapsed to their best, "
        f"{len(labels) - n_hit_labels} label(s) with no known-repeat hit\n"
    )


if __name__ == "__main__":
    main()
