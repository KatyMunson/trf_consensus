#!/usr/bin/env python
"""
cross_motif_comparison.py

Compare every pair of ranked cluster consensus sequences (across all period
bins AND all ranks within a bin) to flag relationships: essentially the
same repeat unit (rotated/on either strand/slightly different TRF-called
period), or one sequence built from tandem copies of a smaller one.

COVERAGE-AWARE TILING, not a single lump alignment. An earlier version
tested "does k copies of the shorter sequence, concatenated, align well
against the longer one" as ONE alignment and used the aggregate identity.
That's a real bug: if the shorter sequence is a small, near-perfect
embedded cassette that only makes up part of the longer sequence, the
aggregate identity over the whole concatenated alignment can still clear
the threshold even when a large fraction of the longer sequence is
completely unexplained — averaging hides exactly the case you most need to
catch (a big fraction of the target being real, distinct content).

Instead, for every pair this greedily tiles the target with non-overlapping
copies of the shorter sequence (rotation- and strand-aware via edlib infix
search on a doubled target, same trick as before), masking each match
before searching for the next, and reports:
  - coverage: fraction of the longer sequence's own length actually
    explained by matched copies (the key new safeguard)
  - mean_copy_identity: average identity of the copies found
  - n_copies_found: how many non-overlapping copies were tiled
A relationship's initial pass/fail comes from coverage + mean identity
thresholds, but that alone isn't sufficient — giving a pair up to
--max-copies independent chances to find a match inflates the false-
positive rate relative to what a flat identity threshold implies (repeated
independent trials, each with some non-zero chance of clearing threshold by
composition alone, multiply up). So every pair that passes the flat
thresholds gets a SECOND check: the same tiling procedure is re-run
--n-shuffles times against a mononucleotide-shuffled version of the unit
(same length, same base composition, order destroyed), and the real
coverage must beat the shuffled null — by default, ALL shuffled trials must
fail to reach --min-coverage. This is a real randomization significance
test, not just a stricter identity cutoff. n_copies_found==1 is reported as
"similar" (same underlying repeat, roughly the same length); n_copies_found
>=2 is reported as "multiple_of" (the longer one is a higher-order/tandem-
duplicate build from the shorter one).
"""
import argparse
import itertools
import random
import sys
from collections import defaultdict

import edlib

COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
    "TGCAYRSWMKVHDBNtgcayrswmkvhdbn",
)


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def read_multi_fasta(paths):
    """Read all records across one or more (possibly multi-record) FASTA
    files. Label is each record's header up to the first space."""
    entries = []
    for path in paths:
        header, seq = None, []
        with open(path) as f:
            for line in f:
                line = line.rstrip("\n")
                if not line:
                    continue
                if line.startswith(">"):
                    if header is not None and seq:
                        entries.append((header, "".join(seq).upper()))
                    header = line[1:].split()[0] if line[1:].split() else line[1:].strip()
                    seq = []
                else:
                    seq.append(line.strip())
            if header is not None and seq:
                entries.append((header, "".join(seq).upper()))
    return entries


def tile_coverage(unit, target, min_identity, max_copies):
    """Greedily tile `target` with non-overlapping copies of `unit`,
    rotation-tolerant (doubled target) and strand-aware, masking each match
    before searching for the next. Tries both strands, keeps whichever
    achieves higher coverage. Returns (coverage_fraction, mean_identity,
    n_copies_found)."""
    best = (0.0, 0.0, 0)
    n = len(target)
    if n == 0 or len(unit) == 0:
        return best
    for strand_target in (target, revcomp(target)):
        working = list(strand_target + strand_target)
        total_matched = 0
        identities = []
        for _ in range(max_copies):
            current = "".join(working)
            result = edlib.align(unit, current, mode="HW", task="path")
            if result["editDistance"] < 0:
                break
            identity = 1.0 - (result["editDistance"] / max(len(unit), 1))
            if identity < min_identity:
                break
            start, end = result["locations"][0]
            span = end - start + 1
            total_matched += span
            identities.append(identity)
            for i in range(start, min(end + 1, len(working))):
                working[i] = "N"
                # `working` is the target doubled to support rotation, so
                # every physical base appears at both i and its mirror i+/-n.
                # Mask both, or an already-counted copy can be "found" again
                # via its mirror image on a later iteration and inflate
                # total_matched / n_copies_found for the same real bases.
                mirror = i - n if i >= n else i + n
                if 0 <= mirror < len(working):
                    working[mirror] = "N"
            if total_matched >= n:
                break
        coverage = min(total_matched / n, 1.0)
        mean_identity = sum(identities) / len(identities) if identities else 0.0
        if coverage > best[0]:
            best = (coverage, mean_identity, len(identities))
    return best


def dinucleotide_shuffle(seq, rng, max_tries=100):
    """Shuffle `seq` preserving its dinucleotide (and single-base)
    composition, via the Altschul-Erikson (1985) edge-shuffle algorithm:
    fix each character's LAST outgoing transition to what it was in the
    original sequence (which guarantees, by construction, that every
    character can still reach the sequence's final character by always
    following its own fixed last edge), shuffle only the order of each
    character's other outgoing transitions, then re-walk the graph.

    EXPERIMENTAL: available via --shuffle-mode di but not the default —
    mononucleotide shuffling only controls for single-base composition,
    which can be too permissive a null for compositionally structured
    satellite sequence (short internal runs, CpG depletion, etc). Needs
    validation against the existing mononucleotide null on real data before
    switching the pipeline default.
    """
    chars = list(seq)
    if len(chars) < 3:
        rng.shuffle(chars)
        return "".join(chars)

    last_char = chars[-1]
    edges = defaultdict(list)
    for i in range(len(chars) - 1):
        edges[chars[i]].append(chars[i + 1])

    def reachable_via_last_edge(shuffled):
        for c, nxts in shuffled.items():
            if c == last_char or not nxts:
                continue
            cur, seen = c, set()
            while cur != last_char:
                if cur in seen or cur not in shuffled or not shuffled[cur]:
                    return False
                seen.add(cur)
                cur = shuffled[cur][-1]
        return True

    shuffled = None
    for _ in range(max_tries):
        candidate = {}
        for c, nxts in edges.items():
            if len(nxts) > 1:
                head = nxts[:-1]
                rng.shuffle(head)
                candidate[c] = head + [nxts[-1]]
            else:
                candidate[c] = list(nxts)
        if reachable_via_last_edge(candidate):
            shuffled = candidate
            break

    if shuffled is None:
        # The fixed-last-edge subgraph is connected to last_char by
        # construction, so this should never actually trigger; fall back to
        # a mononucleotide shuffle rather than fail the whole comparison.
        rng.shuffle(chars)
        return "".join(chars)

    pos = {c: 0 for c in shuffled}
    out = [chars[0]]
    cur = chars[0]
    for _ in range(len(chars) - 1):
        nxt = shuffled[cur][pos[cur]]
        pos[cur] += 1
        out.append(nxt)
        cur = nxt
    return "".join(out)


def null_pass_count(unit, target, min_identity, max_copies, min_coverage, n_shuffles, rng,
                     shuffle_mode="mono"):
    """How many of n_shuffles shuffled versions of `unit` achieve
    >=min_coverage against the same target via the same tiling procedure.
    Real matches should score far better than any shuffle. shuffle_mode
    'mono' (default) preserves single-base composition only; 'di' preserves
    dinucleotide composition too (see dinucleotide_shuffle)."""
    n_pass = 0
    unit_list = list(unit)
    for _ in range(n_shuffles):
        if shuffle_mode == "di":
            shuffled = dinucleotide_shuffle(unit, rng)
        else:
            rng.shuffle(unit_list)
            shuffled = "".join(unit_list)
        cov, _, _ = tile_coverage(shuffled, target, min_identity, max_copies)
        if cov >= min_coverage:
            n_pass += 1
    return n_pass


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fastas", nargs="+", required=True,
                     help="One or more FASTA files (each may hold multiple records); "
                          "each record's own header (up to the first space) is used as its label")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--min-identity", type=float, default=0.70,
                     help="Per-copy identity threshold")
    ap.add_argument("--min-coverage", type=float, default=0.85,
                     help="Fraction of the longer sequence's length that must be explained by "
                          "matched copies of the shorter one for a relationship to be flagged")
    ap.add_argument("--max-copies", type=int, default=6,
                     help="Largest number of tiled copies to search for per pair")
    ap.add_argument("--n-shuffles", type=int, default=20,
                     help="Mononucleotide-shuffled null trials per candidate-passing pair")
    ap.add_argument("--max-null-passes", type=int, default=0,
                     help="A pair is only confirmed if at most this many of --n-shuffles shuffled "
                          "trials also clear --min-coverage. Default 0 = ALL shuffles must fail.")
    ap.add_argument("--shuffle-seed", type=int, default=1)
    ap.add_argument("--shuffle-mode", choices=["mono", "di"], default="mono",
                     help="EXPERIMENTAL, not yet wired into the Snakemake config — needs "
                          "validation before switching the pipeline default. Null-model shuffle "
                          "for the randomization test: 'mono' (default) preserves single-base "
                          "composition only; 'di' is an Altschul-Erikson dinucleotide-preserving "
                          "shuffle, a stricter null more appropriate for compositionally "
                          "structured satellite sequence.")
    args = ap.parse_args()

    entries = read_multi_fasta(args.fastas)
    if not entries:
        sys.stderr.write("[cross_motif_comparison] no sequences found in input\n")

    rng = random.Random(args.shuffle_seed)

    columns = [
        "label_A", "label_B", "len_A", "len_B", "length_ratio",
        "coverage", "mean_copy_identity", "n_copies_found",
        "null_passes", "flag_similar", "flag_multiple_of",
    ]
    rows = []
    n_candidate_pairs = 0

    for (label_a, seq_a), (label_b, seq_b) in itertools.combinations(entries, 2):
        # Order so A is the shorter (or equal) sequence — A is always the
        # tiling unit, B is always the target being tiled/explained.
        if len(seq_a) > len(seq_b):
            label_a, seq_a, label_b, seq_b = label_b, seq_b, label_a, seq_a
        len_a, len_b = len(seq_a), len(seq_b)
        ratio = len_b / len_a if len_a else float("nan")

        max_copies = min(args.max_copies, int(ratio) + 2) if len_a else 1
        coverage, mean_identity, n_copies = tile_coverage(
            seq_a, seq_b, args.min_identity, max_copies
        )

        candidate = (coverage >= args.min_coverage) and (mean_identity >= args.min_identity) and n_copies >= 1
        null_passes = "NA"
        flag_similar = flag_multiple = False

        if candidate:
            n_candidate_pairs += 1
            null_passes = null_pass_count(
                seq_a, seq_b, args.min_identity, max_copies, args.min_coverage,
                args.n_shuffles, rng, shuffle_mode=args.shuffle_mode
            )
            confirmed = null_passes <= args.max_null_passes
            flag_similar = confirmed and n_copies == 1
            flag_multiple = confirmed and n_copies >= 2

        rows.append({
            "label_A": label_a, "label_B": label_b,
            "len_A": len_a, "len_B": len_b,
            "length_ratio": f"{ratio:.3f}",
            "coverage": f"{coverage:.4f}",
            "mean_copy_identity": f"{mean_identity:.4f}" if n_copies else "NA",
            "n_copies_found": n_copies,
            "null_passes": null_passes,
            "flag_similar": flag_similar,
            "flag_multiple_of": flag_multiple,
        })

    with open(args.out_tsv, "w") as out:
        out.write("\t".join(columns) + "\n")
        for row in rows:
            out.write("\t".join(str(row[c]) for c in columns) + "\n")

    n_flagged = sum(1 for r in rows if r["flag_similar"] or r["flag_multiple_of"])
    sys.stderr.write(
        f"[cross_motif_comparison] {len(entries)} sequences, {len(rows)} pairs compared, "
        f"{n_candidate_pairs} passed the flat threshold, {n_flagged} confirmed by "
        f"randomization test (max {args.max_null_passes}/{args.n_shuffles} shuffles allowed to pass)\n"
    )


if __name__ == "__main__":
    main()
