#!/usr/bin/env python
"""
rank_family_clusters.py

For one period bin (01_raw_consensus.fasta), repeatedly: pick a reference
(median-length sequence among whatever remains), anchor-match everyone else
against it (same rotation/strand-aware edlib logic as the old
canonicalize_rotation.py), peel off what matches as one cluster, and repeat
on the leftovers. Stops once the remaining pool (or the newly found
cluster) is smaller than --min-cluster-size, or --max-rounds is hit.
Clusters are then sorted by size and reported as rank 1 (most loci), rank
2, etc. Sequences that never join a cluster of the minimum size are
dropped (not written anywhere) — for centromere/large-satellite hunting the
signal-bearing clusters are what matter, not the noise tail.

Why iterative peeling instead of full all-vs-all clustering: all-vs-all is
O(N^2) pairwise alignments — for a bin with a few thousand loci that's
millions of pairs, hours of compute on one core. This is O(rounds x N):
each round is one reference vs. the remaining pool (an O(N) anchor-match
pass), and the number of rounds is just however many distinct coherent
families actually exist in the bin (usually small). It also makes the old
"primary/secondary swap" unnecessary — nothing is labeled "primary" until
every cluster has been found and sizes compared.

For each surviving cluster, this calls mafft for the MSA and reimplements
the same majority-rule consensus caller as the old consensus_from_alignment.py
(duplicated inline to keep this script self-contained — the cluster count
per bin isn't known ahead of time, so chaining separate Snakemake jobs per
cluster isn't practical here).

Outputs, all inside --out-dir:
  rank01_n<N>.fasta, rank02_n<N>.fasta, ... — one consensus per cluster,
    header ">{motif_name}_rank{rank}_n{N} n_input_sequences=N consensus_length=L"
  --out-summary: one row per cluster, sorted by rank
"""
import argparse
import os
import random
import subprocess
import sys
import tempfile
from collections import Counter

import edlib

COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
    "TGCAYRSWMKVHDBNtgcayrswmkvhdbn",
)
ANCHOR_LEN = 21
ANCHOR_START_OFFSETS = [0.5, 0.25, 0.75, 0.1, 0.9]


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


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


def best_anchor_hit(anchor, target_doubled, min_identity):
    result = edlib.align(anchor, target_doubled, mode="HW", task="locations")
    if result["editDistance"] < 0:
        return None, 0.0
    identity = 1.0 - (result["editDistance"] / len(anchor))
    if identity < min_identity:
        return None, 0.0
    return result["locations"][0][0], identity


def rotate_to_offset(seq, hit_start, ref_anchor_offset):
    n = len(seq)
    shift = (hit_start - ref_anchor_offset) % n
    doubled = seq + seq
    return doubled[shift:shift + n]


def anchor_is_unique(anchor, ref_seq, own_start, min_identity, ambiguity_margin=0.05):
    """Defensive check for --require-unique-anchor: True unless `anchor`
    (taken from ref_seq at own_start) has another equally-good match
    elsewhere in ref_seq. Guards against the reference sitting in an
    internally repeated region (e.g. higher-order-repeat structure), where
    anchor-based rotation could lock different cluster members onto
    different internal copies instead of the same homologous position.
    edlib's "locations" task only reports optimal (tied-best) hits, so this
    catches exact/near-exact internal duplication, not every possible
    ambiguity — it's a heuristic, not a complete check.
    """
    doubled = ref_seq + ref_seq
    result = edlib.align(anchor, doubled, mode="HW", task="locations")
    if result["editDistance"] < 0:
        return True
    identity = 1.0 - (result["editDistance"] / len(anchor))
    if identity < min_identity - ambiguity_margin:
        return True
    n = len(ref_seq)
    other_hits = [loc for loc in result["locations"] if loc[0] % n != own_start]
    return len(other_hits) == 0


def peel_one_cluster(records, min_identity, anchor_len=ANCHOR_LEN,
                      anchor_offsets=ANCHOR_START_OFFSETS, require_unique_anchor=False):
    """One round: pick a reference (median length), anchor-match everyone
    else against it. Returns (cluster_records_rotated, leftover_records)."""
    lengths = sorted(len(s) for _, s in records)
    median_len = lengths[len(lengths) // 2]
    ref_id, ref_seq = min(records, key=lambda r: abs(len(r[1]) - median_len))

    cluster = [(ref_id, ref_seq)]
    leftover = []

    for rec_id, seq in records:
        if rec_id == ref_id:
            continue
        best = None
        for offset_frac in anchor_offsets:
            anchor_start = int(len(ref_seq) * offset_frac)
            anchor = ref_seq[anchor_start:anchor_start + anchor_len]
            if len(anchor) < anchor_len:
                continue
            if require_unique_anchor and not anchor_is_unique(anchor, ref_seq, anchor_start, min_identity):
                continue
            for strand, cand in (("+", seq), ("-", revcomp(seq))):
                doubled = cand + cand
                hit_start, identity = best_anchor_hit(anchor, doubled, min_identity)
                if hit_start is None:
                    continue
                rotated = rotate_to_offset(cand, hit_start, anchor_start)
                if best is None or identity > best[0]:
                    best = (identity, rotated)
            if best is not None and best[0] >= 0.95:
                break
        if best is None:
            leftover.append((rec_id, seq))
        else:
            cluster.append((rec_id, best[1]))

    return cluster, leftover


def run_mafft(records, mafft_bin):
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        for rec_id, seq in records:
            f.write(f">{rec_id}\n{seq}\n")
        in_path = f.name
    try:
        result = subprocess.run([mafft_bin, "--auto", in_path], capture_output=True, text=True)
    finally:
        os.unlink(in_path)
    if result.returncode != 0:
        sys.stderr.write(f"WARNING: mafft failed: {result.stderr[:500]}\n")
        return None
    return result.stdout


def read_fasta_text(text):
    ids, seqs = [], []
    seq = []
    for line in text.splitlines():
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


def majority_consensus(aligned_fasta_text, max_gap_fraction, rng=None):
    """Majority-rule consensus per alignment column. Ties (no single most-
    common base) are broken with `rng` rather than always favoring
    A>C>G>T, which would otherwise introduce a small systematic base-
    composition bias into every consensus with any tied columns."""
    records = read_fasta_text(aligned_fasta_text)
    seqs = [s.upper() for _, s in records]
    if not seqs:
        return ""
    if rng is None:
        rng = random.Random(0)
    aln_len = len(seqs[0])
    n = len(seqs)
    bases = []
    for col in range(aln_len):
        column = [s[col] for s in seqs]
        if column.count("-") / n >= max_gap_fraction:
            continue
        counts = Counter(b for b in column if b in "ACGT")
        if not counts:
            continue
        top = counts.most_common()
        max_count = top[0][1]
        tied = sorted(b for b, c in top if c == max_count)
        bases.append(tied[0] if len(tied) == 1 else rng.choice(tied))
    return "".join(bases)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in-fasta", required=True, help="01_raw_consensus.fasta for one motif bin")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--motif-name", required=True)
    ap.add_argument("--min-identity", type=float, default=0.80)
    ap.add_argument("--max-gap-fraction", type=float, default=0.5)
    ap.add_argument("--min-cluster-size", type=int, default=3,
                     help="Stop peeling once the remaining pool (or a newly found cluster) is smaller than this")
    ap.add_argument("--max-rounds", type=int, default=15)
    ap.add_argument("--mafft-bin", default="mafft")
    ap.add_argument("--consensus-tie-seed", type=int, default=1,
                     help="Seed for random tie-breaking of majority-consensus columns with no "
                          "single most-common base (previously always broke ties toward "
                          "A>C>G>T, a small systematic composition bias)")
    ap.add_argument("--anchor-length", type=int, default=ANCHOR_LEN,
                     help="EXPERIMENTAL knob, not yet wired into the Snakemake config: anchor "
                          "length used for rotation/strand matching. Default reproduces the "
                          "prior hardcoded behavior.")
    ap.add_argument("--anchor-offsets", default=",".join(str(x) for x in ANCHOR_START_OFFSETS),
                     help="EXPERIMENTAL knob, not yet wired into the Snakemake config: comma-"
                          "separated fractional offsets (0-1) into the reference sequence to try "
                          "anchors at, in order. Default reproduces the prior hardcoded behavior.")
    ap.add_argument("--require-unique-anchor", action="store_true",
                     help="EXPERIMENTAL, off by default and not yet wired into the Snakemake "
                          "config — needs validation on real data before becoming the default. "
                          "Rejects an anchor that has another equally-good match elsewhere in "
                          "the same reference before using it for rotation, guarding against "
                          "locking onto the wrong copy in an internally repeated (e.g. higher-"
                          "order-repeat) reference.")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    anchor_offsets = [float(x) for x in args.anchor_offsets.split(",")]
    consensus_rng = random.Random(args.consensus_tie_seed)

    records = read_fasta(args.in_fasta)
    n_total = len(records)
    sys.stderr.write(f"[rank_family_clusters] {args.motif_name}: {n_total} input sequences\n")

    clusters = []
    remaining = records
    round_n = 0
    n_dropped_small = 0
    while remaining and len(remaining) >= args.min_cluster_size and round_n < args.max_rounds:
        round_n += 1
        cluster, remaining = peel_one_cluster(
            remaining, args.min_identity,
            anchor_len=args.anchor_length,
            anchor_offsets=anchor_offsets,
            require_unique_anchor=args.require_unique_anchor,
        )
        if len(cluster) < args.min_cluster_size:
            # This round's reference didn't anchor a coherent-enough
            # cluster. Drop just these sequences and keep trying with the
            # rest of the pool — a single unlucky median-length reference
            # shouldn't hide a real, differently-sized family still sitting
            # in `remaining` (previously this `break`d and discarded the
            # entire remaining pool on the first bad reference).
            n_dropped_small += len(cluster)
            sys.stderr.write(
                f"[rank_family_clusters] {args.motif_name} round {round_n}: reference only "
                f"anchored {len(cluster)} (<{args.min_cluster_size}), dropping and continuing "
                f"with {len(remaining)} remaining\n"
            )
            continue
        clusters.append(cluster)
        sys.stderr.write(
            f"[rank_family_clusters] {args.motif_name} round {round_n}: peeled off "
            f"{len(cluster)}, {len(remaining)} remain\n"
        )

    n_unclustered = len(remaining) + n_dropped_small
    clusters.sort(key=len, reverse=True)  # rank 1 = largest = most support

    summary_rows = []
    n_consensus_failed = 0
    for rank, cluster in enumerate(clusters, 1):
        aligned_text = run_mafft(cluster, args.mafft_bin)
        consensus = majority_consensus(aligned_text, args.max_gap_fraction, consensus_rng) if aligned_text else ""

        if not consensus:
            n_consensus_failed += 1
            sys.stderr.write(
                f"ERROR: {args.motif_name} rank{rank:02d} (n={len(cluster)}) produced an empty "
                f"consensus (mafft failure or all-gap alignment) — dropping this cluster rather "
                f"than emitting a phantom record with no sequence\n"
            )
            continue

        label = f"{args.motif_name}_rank{rank}_n{len(cluster)}"
        out_path = os.path.join(args.out_dir, f"rank{rank:02d}_n{len(cluster)}.fasta")
        with open(out_path, "w") as f:
            f.write(f">{label} n_input_sequences={len(cluster)} consensus_length={len(consensus)}\n")
            f.write(consensus + "\n")

        gc = (consensus.count("G") + consensus.count("C")) / len(consensus) if consensus else 0.0
        summary_rows.append({
            "label": label, "rank": rank, "n_input_sequences": len(cluster),
            "pct_of_bin": f"{100.0 * len(cluster) / n_total:.2f}" if n_total else "NA",
            "consensus_length": len(consensus), "gc_content": f"{gc:.4f}", "fasta": out_path,
        })

    with open(args.out_summary, "w") as out:
        out.write("label\trank\tn_input_sequences\tpct_of_bin\tconsensus_length\tgc_content\tfasta\n")
        for r in summary_rows:
            out.write(f"{r['label']}\t{r['rank']}\t{r['n_input_sequences']}\t{r['pct_of_bin']}\t"
                      f"{r['consensus_length']}\t{r['gc_content']}\t{r['fasta']}\n")

    sys.stderr.write(
        f"[rank_family_clusters] {args.motif_name}: {len(clusters) - n_consensus_failed} clusters "
        f"written in {round_n} rounds, {n_unclustered}/{n_total} sequences dropped as unclustered, "
        f"{n_consensus_failed} cluster(s) dropped for empty/failed consensus\n"
    )


if __name__ == "__main__":
    main()
