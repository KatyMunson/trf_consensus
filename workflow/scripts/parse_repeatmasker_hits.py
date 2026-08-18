#!/usr/bin/env python
"""
parse_repeatmasker_hits.py

Turn RepeatMasker's .out file (from the known-repeat screening step, run
against our own discovered consensus motifs) into a clean per-label TSV.

Every label in --query-fasta is represented, whether or not RepeatMasker
found anything for it — RepeatMasker's .out file only ever lists sequences
with at least one hit, so the full label universe has to come from the
query FASTA itself, not from the .out file. A missing .out file
(RepeatMasker can decline to create one at all when nothing anywhere in the
query matches anything in the library) or a header-only one is treated
exactly like zero hits found, not an error: for satellite DNA in a
non-model species this is a completely plausible, even expected, real
result. A label with no hit gets exactly one row (has_known_hit=False, rest
NA). A label WITH hits gets one row per DISTINCT matched repeat name — see
"Multiple hits" below — so this table is not guaranteed one-row-per-label.

Parsing approach: RepeatMasker's .out format has a few header/blank lines
before data rows begin, and the exact header formatting has drifted
slightly across versions. Rather than hardcoding a line count, a line is
treated as a data row if, after stripping and splitting on whitespace, its
first token parses as an integer (the SW score column) — the same heuristic
already used in extract_by_period.py for skipping TRF .dat header lines.

Column layout (whitespace-delimited). The three "position in repeat"
columns swap order depending on strand: for "+" they're
(begin, end, (left)); for "C" (minus strand) RepeatMasker reports them as
((left), end, begin) instead — same physical meaning, different column
order, so this is parsed strand-aware rather than by fixed position:
  [0] SW score            [8]  strand (+ or C)
  [1] % divergence        [9]  matching repeat name
  [2] % deleted            [10] repeat class/family
  [3] % inserted           [11] position in repeat: begin, or (left) if strand==C
  [4] query sequence name  [12] position in repeat: end
  [5] query begin          [13] position in repeat: (left), or begin if strand==C
  [6] query end             [14] ID
  [7] query (left), parenthesized
Strand is normalized from RepeatMasker's "C" to "-" for consistency with
this pipeline's own strand convention elsewhere ("+"/"-").

Overlap columns (derived, not RepeatMasker's own output):
  - pct_query_covered: fraction of THIS query's true (un-doubled) length
    spanned by this hit's query_start/query_end. If --doubled was used,
    the raw span is capped at the true (pre-doubling) motif length before
    dividing, the same way cross_motif_comparison.py's tile_coverage()
    avoids double-counting via the doubling trick — otherwise a hit that
    wraps into the second copy could nonsensically exceed 100%.
  - pct_known_repeat_covered: fraction of the MATCHED REPEAT's own total
    model/consensus length spanned by this hit, derived from RepeatMasker's
    "position in repeat" begin/end/left columns (end - begin + 1 explained,
    out of a total model length of end + left).
  - reciprocal_overlap: min(pct_query_covered, pct_known_repeat_covered) —
    the conservative combined metric; both sequences must be substantially
    explained for this to be high. Guards against e.g. a huge query that
    happens to embed a tiny fragment of a known repeat (high
    pct_known_repeat_covered, low pct_query_covered) reading as a strong
    hit when only a small piece is actually explained, or the reverse.

Multiple hits per label: queries are doubled (seq+seq) when --doubled is
set (matching known_repeat_screen.double_sequences), so up to 2 raw hit
rows can come back per true label for what's really the SAME underlying
match — one from each copy in the doubled query. These collapse to one
output row by keeping the hit with the highest SW score PER (label,
matched repeat name) — so a genuinely distinct second match, against a
DIFFERENT known repeat family, is kept as its own separate row rather than
being discarded in favor of whichever hit happened to score higher. This
means a label can appear multiple times in this table; if you need exactly
one row per label (e.g. for a simple join against summary_table.tsv), sort
by reciprocal_overlap descending and keep the first row per label.

CAVEAT, not something this script tries to correct for: if a real match
happens to span the artificial junction introduced by doubling, RepeatMasker
can report an inflated or split score right at that boundary. This is a
known, accepted approximation of the doubling trick. If a top hit looks
borderline or surprising, it's worth a manual look at the un-doubled
first-half coordinates (query_start/query_end here are reported exactly as
RepeatMasker output them, i.e. in the doubled query's coordinate space when
--doubled was used — only the derived pct_query_covered column corrects
for doubling).
"""
import argparse
import os
import sys

STRAND_MAP = {"+": "+", "C": "-"}


def read_fasta_labels_and_lengths(path):
    """Header labels in file order (the complete label universe the output
    table must cover) plus each record's sequence length. Assumes one
    sequence line per record (this pipeline's own FASTA convention, see
    prepare_known_repeat_query.py), but accumulates across lines regardless
    in case that ever changes."""
    labels = []
    lengths = {}
    label = None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                label = line[1:].strip().split()[0]
                labels.append(label)
                lengths[label] = 0
            elif label is not None:
                lengths[label] += len(line.strip())
    return labels, lengths


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
            if len(fields) < 14 or not fields[0].isdigit():
                continue
            try:
                strand = STRAND_MAP.get(fields[8], fields[8])
                if fields[8] == "C":
                    rep_left, rep_end, rep_begin = fields[11], fields[12], fields[13]
                else:
                    rep_begin, rep_end, rep_left = fields[11], fields[12], fields[13]
                yield {
                    "sw_score": int(fields[0]),
                    "pct_divergence": fields[1],
                    "pct_deleted": fields[2],
                    "pct_inserted": fields[3],
                    "query_label": fields[4],
                    "query_start": int(fields[5]),
                    "query_end": int(fields[6]),
                    "strand": strand,
                    "repeat_name": fields[9],
                    "repeat_class_family": fields[10],
                    "repeat_begin": int(rep_begin),
                    "repeat_end": int(rep_end),
                    "repeat_left": int(rep_left.strip("()")),
                }
            except (IndexError, ValueError):
                continue


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rm-out", required=True, help="RepeatMasker .out file (may not exist)")
    ap.add_argument("--query-fasta", required=True,
                     help="The query FASTA that was screened (prepare_known_repeat_query.py's "
                          "output) — source of the full label universe and true motif lengths, "
                          "since RepeatMasker's .out only lists sequences with a hit")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--doubled", action="store_true",
                     help="Matches known_repeat_screen.double_sequences — used to recover each "
                          "motif's true (pre-doubling) length for pct_query_covered")
    args = ap.parse_args()

    labels, lengths = read_fasta_labels_and_lengths(args.query_fasta)
    true_length = {
        label: (length // 2 if args.doubled else length)
        for label, length in lengths.items()
    }

    # best hit per (label, matched repeat name) — collapses the doubling
    # artifact (same real match, found once per copy) without discarding a
    # genuinely distinct match against a different known repeat family
    best_hit = {}
    n_rows = 0
    for row in parse_rm_out(args.rm_out):
        n_rows += 1
        key = (row["query_label"], row["repeat_name"])
        if key not in best_hit or row["sw_score"] > best_hit[key]["sw_score"]:
            best_hit[key] = row

    hits_by_label = {}
    for (label, _repeat_name), hit in best_hit.items():
        hits_by_label.setdefault(label, []).append(hit)

    columns = ["label", "has_known_hit", "repeat_name", "repeat_class_family", "sw_score",
               "pct_divergence", "pct_deleted", "pct_inserted", "query_start", "query_end", "strand",
               "pct_query_covered", "pct_known_repeat_covered", "reciprocal_overlap"]

    rows_out = []
    for label in labels:
        hits = sorted(hits_by_label.get(label, []), key=lambda h: h["sw_score"], reverse=True)
        if not hits:
            row = {c: "NA" for c in columns}
            row["label"] = label
            row["has_known_hit"] = False
            rows_out.append(row)
            continue

        motif_len = true_length.get(label, 0)
        for hit in hits:
            query_span = hit["query_end"] - hit["query_start"] + 1
            pct_query = min(query_span, motif_len) / motif_len if motif_len else 0.0

            repeat_model_len = hit["repeat_end"] + hit["repeat_left"]
            repeat_span = hit["repeat_end"] - hit["repeat_begin"] + 1
            pct_repeat = repeat_span / repeat_model_len if repeat_model_len else 0.0

            rows_out.append({
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
                "pct_query_covered": f"{pct_query:.4f}",
                "pct_known_repeat_covered": f"{pct_repeat:.4f}",
                "reciprocal_overlap": f"{min(pct_query, pct_repeat):.4f}",
            })

    with open(args.out_tsv, "w") as out:
        out.write("\t".join(columns) + "\n")
        for row in rows_out:
            out.write("\t".join(str(row[c]) for c in columns) + "\n")

    n_hit_labels = len(hits_by_label)
    n_distinct_hits = len(best_hit)
    sys.stderr.write(
        f"[parse_repeatmasker_hits] {len(labels)} query labels, {n_rows} raw RepeatMasker hit "
        f"row(s) parsed, {n_hit_labels} label(s) with >=1 hit ({n_distinct_hits} distinct "
        f"label/repeat matches, {len(rows_out)} output rows), "
        f"{len(labels) - n_hit_labels} label(s) with no known-repeat hit\n"
    )


if __name__ == "__main__":
    main()
