# TRF Motif Consensus-of-Consensus Pipeline

Fully automatic, no manually-curated candidate list: scans a raw TRF `.dat`
file for candidate period bins, builds a ranked set of consensus sequences
per bin via iterative clustering, resolves redundancy across every cluster
found (including tandem-duplicate relationships between different bins),
and outputs a single RepeatMasker-ready custom library plus a QC summary
table.

## Pipeline

```
scan_dat_candidates (checkpoint)
        |
        v
extract_by_period  (per bin)
        |
        v
rank_family_clusters  (per bin: iterative clustering + per-cluster consensus)
        |
        v
build_all_clusters_table  +  combine_cluster_fastas
        |
        v
cross_cluster_comparison
        |
        v
resolve_redundancy
        |
        v
results/repeatmasker_custom_lib.fasta + results/summary_table.tsv
```

1. **`scan_dat_candidates`** scans the whole `.dat` file directly and
   nominates every period bin meeting the `candidate_scan` filters (see
   below). This is a Snakemake **checkpoint** — the bin list isn't known
   until this actually runs, so everything downstream is generated
   dynamically from its output, not from a static manifest file. There's no
   manual pause or curation step: every passing bin proceeds automatically.
2. **`extract_by_period`** pulls each bin's loci out of the `.dat` file
   (unchanged from earlier versions of this pipeline).
3. **`rank_family_clusters`** replaces the old two-round primary/secondary
   system entirely. For each bin, it repeatedly: picks a reference sequence
   (median length among whatever remains), anchor-matches everyone else
   against it (rotation- and strand-aware, via `edlib`), peels off what
   matches as one cluster, and repeats on the leftovers — until the
   remaining pool is smaller than `min_cluster_size` or `max_rounds` is
   hit. All discovered clusters in a bin are then sorted by size: rank 1 =
   most support, rank 2 = next, etc. Each cluster gets its own MAFFT
   alignment and majority-rule consensus. **Sequences that never join a
   cluster of the minimum size are dropped, not reported** — this pipeline
   is tuned for finding high-signal candidates (centromeric/large satellite
   arrays), not for characterizing background noise.

   This is deliberately **not** full all-vs-all pairwise clustering.
   All-vs-all is O(N²) alignments — a bin with a few thousand loci is
   millions of pairs, hours of compute on one core. Iterative peeling is
   O(rounds × N): each round is one reference vs. the remaining pool (an
   O(N) pass), and the round count is just however many distinct coherent
   families actually exist in the bin (usually small). It also makes the
   old primary/secondary "swap" logic unnecessary — nothing is labeled
   primary until every cluster in every bin has been found and sizes
   compared globally, in the redundancy-resolution step below.
4. **`build_all_clusters_table`** + **`combine_cluster_fastas`** aggregate
   every bin's ranked clusters into one master table and one combined FASTA.
5. **`cross_cluster_comparison`** compares every cluster found — all bins,
   all ranks, including different ranks within the same bin — using
   **coverage-aware tiling with a randomization significance test**, not a
   single lump alignment. For each pair (shorter sequence A, longer
   sequence B), it greedily tiles B with non-overlapping copies of A
   (rotation- and strand-aware, masking each match before searching for the
   next) and requires both a per-copy identity threshold AND that the tiled
   copies explain most of B's length (`min_coverage`, default 0.85) before
   even considering a relationship. This replaced an earlier version that
   concatenated exactly *k* copies of A and aligned the whole thing as one
   block — averaging identity over that whole alignment let a small,
   near-perfect *embedded cassette* match drag up the aggregate score even
   when most of B was genuinely unrelated content, silently discarding real
   distinct motifs that happen to share a common building-block subsequence
   with something more abundant. Every pair that clears the coverage/identity
   bar then gets a **second, statistical check**: the same tiling procedure
   is re-run against `n_shuffles` (default 20) mononucleotide-shuffled
   versions of A, and the real match must beat all (or `max_null_passes`)
   of those shuffled attempts — this catches cases where repeatedly giving a
   pair several independent tiling attempts (up to `max_copies`) would
   otherwise inflate the false-positive rate beyond what the raw identity
   threshold implies, especially for short or compositionally simple
   sequences.
6. **`resolve_redundancy`** processes clusters in descending support order
   (highest `n_input_sequences` first). Each cluster is checked only
   against clusters **already confirmed as winners** — never against
   another loser, never chained through an intermediate. If it's flagged
   against an existing winner, it's marked redundant with that winner
   specifically (a direct, individually-verified relationship); otherwise
   it becomes a new winner itself. **This deliberately avoids transitive
   union-find** (A~B flagged + B~C flagged does NOT imply A~C) — plain
   union-find over flagged pairs was tried first and, on real data,
   collapsed dozens of genuinely unrelated candidates spanning an 11x
   period range into one meaningless "family" simply because each
   consecutive pair along a chain happened to be individually flagged.
   Nothing is dropped — losers are kept in both the summary table and the
   final RepeatMasker library, marked `is_redundant=True` / demoted to
   `#Satellite/redundant_with_<winner>` in the library's classification
   field, so they're still visible (and RepeatMasker-usable) but clearly
   lower priority.

## Sanity-checking a mega-group

If one cluster ends up "winning" a large redundancy group (many clusters
spanning a wide range of target periods all marked redundant with it),
don't assume that's necessarily wrong — it can be a genuine finding: a
small, highly abundant repeat unit that's been incorporated as a building
block into many otherwise-distinct larger motifs across the genome (a real
pattern for old, dominant satellite families). Before trusting it, check:
`coverage` and `n_copies_found` in `cross_cluster_comparison.tsv` tell you
whether a given partner is *fully* explained by the winner (safe to treat
as redundant) or only *partially* (`coverage` well below 1.0 — that
partner likely has real, distinct content of its own and deserves to be
its own candidate, not demoted). `null_passes` tells you how many of the
20 randomization trials also passed by chance — anything above 0 is worth
a second look.

## Setup

Edit `config/config.yaml`:
- `trf_dat`: path to your `.dat` file
- `trf_dat_default_seqname`: only needed if your `.dat` has no header lines at all
- `repeatmasker_classification`: default `"Satellite"` — edit to match how you
  want these classified before loading into RepeatMasker/FamDB
- **Verify the MAFFT module name** — `rank_family_clusters` requests
  `mafft/7.487` via `envmodules:`. Run `module avail mafft` on liger and
  update the Snakefile if the version differs.
- `known_repeat_screen.species`: the RepeatMasker `-species` value to
  screen discovered motifs against (see "Known-repeat screening" below).

## Run

```bash
snakemake -s Snakefile --configfile config/config.yaml \
    --cores 48 --use-envmodules --use-conda --retries 3
```

**Both `--use-envmodules` and `--use-conda` are required together** — every
rule but one uses envmodules (matching the rest of this pipeline's
cluster-module setup), but `run_repeatmasker_known_screen` deliberately has
no `envmodules:` fallback (it installs its own fresh, pipeline-owned
RepeatMasker via conda rather than reusing whatever's already on the
cluster — see "Known-repeat screening" below). Passing both flags lets
Snakemake use envmodules wherever a rule defines them and fall back to
conda for that one rule; passing `--use-envmodules` alone would run
RepeatMasker with no environment management at all.

(Swap `--use-envmodules` for `--use-conda` to run off the per-rule envs in
`workflow/envs/` instead: `python_base.yaml` for pure-Python rules,
`edlib.yaml` for `cross_cluster_comparison`, `cluster_rank.yaml`
[edlib + mafft together] for `rank_family_clusters`, the only rule that
needs both tools in one process.)

Because `scan_dat_candidates` is a checkpoint, the first `snakemake -n`
dry-run will report "the run involves checkpoint jobs, which will result in
alteration of the DAG of jobs" — that's expected, not an error. The actual
per-bin job count isn't known until the scan completes.

## Config knobs

**`candidate_scan`** — which period bins get processed at all. A bin passes
if EITHER: (distributed) `≥min_blocks` total loci within `period±window`
AND at least one with `copy_number ≥ min_copy_number`; OR (single massive
block) any one locus with `copy_number ≥ min_single_block_copy_number`,
regardless of block count — this catches a genuine large array TRF captured
as one contiguous locus, which the distributed rule alone would reject.
`min_period_length` filters out micro/mini-satellites before either rule is
applied — its default (150) is a human-centric choice (just below alpha-
satellite's 171bp); lower it for other taxa or a permissive discovery pass,
see the comment in `config/configexample.yaml`. `nms_radius` prevents nearby
period values (TRF's period estimate jitters a few bp locus-to-locus for the
same true repeat) from being reported as separate candidates.

**Known gap between the two `candidate_scan` rules:** a genuinely large
satellite array that TRF fragments into many small-to-moderate blocks (e.g.
a HOR array with enough internal degeneracy that no single block reaches
`copy_number ≥ min_single_block_copy_number`, and no individual block
reaches `copy_number ≥ min_copy_number` either) can fall through *both*
rules even though `min_blocks` is satisfied many times over — the
distributed rule requires at least one block to individually clear
`min_copy_number`, not just that the blocks collectively represent a lot of
sequence. If you suspect this is happening for a period you expected to see,
check `total_array_bp` for that period in `candidate_periods_raw.tsv` (every
eligible period's own stats, unsuppressed) — a large `total_array_bp` with
`passes_filters=False` is the signature of this gap, and lowering
`min_copy_number` (rather than `min_blocks`, which is likely already being
met) is the fix.

**`cluster_ranking`** — `min_cluster_size` (default 3) is the floor for a
peeled-off group to count as a real cluster; `max_rounds` (default 15) is a
safety cap; `min_cluster_identity` (default 0.80) is the edlib identity
required to join a cluster; `max_gap_fraction` (default 0.5) is the
alignment-column gap fraction above which a column is dropped from the
consensus.

**`cross_cluster_comparison`** — `min_identity` (default 0.70) is the
per-copy identity required when tiling; `min_coverage` (0.85) is the
fraction of the longer sequence that must be explained by tiled copies;
`max_copies` (6) bounds how many tiled copies are searched for per pair;
`n_shuffles` (20) and `max_null_passes` (0) control the randomization
significance test — raise `max_null_passes` above 0 only if you want a more
permissive (less strict) confirmation criterion.

**`known_repeat_screen`** — `species` is the RepeatMasker `-species` value;
`double_sequences` (default true) controls the rotation-tolerant doubling
trick. See "Known-repeat screening" below for the full rationale.

## Outputs

- `results/candidate_periods.tsv` / `candidate_periods_raw.tsv` /
  `candidate_periods_manifest.tsv` — the scan's output (see
  `scan_dat_candidates.py` docstring for column definitions)
- `results/{motif}/01_raw_consensus.fasta` + `.tsv` — extracted per-locus
  consensus motifs + provenance for that bin
- `results/{motif}/ranked_clusters/rank01_n<N>.fasta`, `rank02_...`, etc. —
  one consensus per cluster found in that bin
- `results/{motif}/ranked_clusters_summary.tsv` — per-bin cluster ranking
- `results/all_clusters.tsv` — every cluster from every bin, one master table
- `results/all_clusters_consensus.fasta` — every cluster's consensus, combined
- `results/cross_cluster_comparison.tsv` — pairwise similarity/multiple-of
  check across every cluster found
- `results/repeatmasker_custom_lib.fasta` — **final output**: every
  cluster's consensus, headers as `>{label}#{classification}` (winners) or
  `>{label}#{classification}/redundant_with_{winner}` (demoted)
- `results/summary_table.tsv` — **final output**: `all_clusters.tsv` plus
  `is_redundant` / `redundant_with` / `group_size`
- `results/known_repeat_hits.tsv` — **final output**: every cluster's
  `label` joined against a Dfam/RepBase known-repeat screen, see
  "Known-repeat screening" below

## Known-repeat screening

A separate, final step screens our own discovered consensus motifs (not
the genome assembly) against a Dfam/RepBase library, via
`RepeatMasker -species`, to check whether any of them match a previously
characterized/named repeat family. This answers a narrower, more direct
question than running RepeatMasker against the whole assembly — "is this
specific candidate we found already known?" — and is much cheaper: ~100-250
short sequences, well under a minute of actual RepeatMasker runtime
regardless of library size, versus screening an entire genome.

This step runs on `results/repeatmasker_custom_lib.fasta` (the final
output — every cluster, winners and redundancy-demoted ones both, since a
demoted cluster individually matching a known family is still useful
information) and is wired into `rule all` alongside the other final
outputs, consistent with this pipeline's fully-automatic design — there's
no separate flag to opt in or out.

**Pipeline**: `prepare_known_repeat_query` (strips our own
`#classification` suffix off each header, since it's not a valid
RepeatMasker query name and `#` has special meaning in RepeatMasker's own
library format) → `run_repeatmasker_known_screen` (RepeatMasker itself) →
`parse_known_repeat_hits` (turns the `.out` file into
`results/known_repeat_hits.tsv`).

**Config** (`known_repeat_screen` in `config/config.yaml`):
- `species` — the RepeatMasker `-species` value. Per FamDB's own docs,
  specifying a species pulls its entire ancestor lineage automatically
  (e.g. Aves, Neognathae, ... all inherited) plus any species-specific
  curated/de novo content filed at that node — specificity is additive
  here, not narrowing, so use the most specific available relative. Before
  relying on a species value, check what's actually populated with
  `famdb.py -i <path-to-famdb> lineage -a "<species>"` (verify the exact
  flag against your installed `famdb.py --help`; this has changed across
  versions — see "One-time setup" below).
- `double_sequences` — concatenates each consensus to itself before
  screening (`seq+seq`), so a database entry that starts at a different
  rotation phase than our consensus can still align end-to-end against a
  contiguous stretch of the query, instead of being missed because the
  real match wraps around the end of a single un-doubled copy. Same trick
  used internally elsewhere in this pipeline (`cross_motif_comparison.py`).

**Reading `known_repeat_hits.tsv`**: every cluster `label` (same join key
as `summary_table.tsv`) is present, but **this is not guaranteed to be one
row per label** — a label with no hit gets exactly one row
(`has_known_hit=False`, rest `NA`), a normal, common, and often *expected*
outcome, not a failure (much of the satellite DNA in a non-model species
genome is genuinely undescribed in existing databases); a label with hits
gets **one row per distinct matched repeat name**, so a motif that matches
two unrelated known families shows up as two rows. If you need exactly one
row per label (e.g. for a simple join), sort by `reciprocal_overlap`
descending and keep the first row per label. `repeat_name` /
`repeat_class_family` identify each match, `sw_score` / `pct_divergence`
describe its quality, and:
- `pct_query_covered` — fraction of *your* motif's true (pre-doubling)
  length this hit explains, capped at 1.0 so a hit wrapping into the
  second copy of a doubled query can't nonsensically exceed 100%.
- `pct_known_repeat_covered` — fraction of the *matched repeat family's
  own* model/consensus length this hit explains, from RepeatMasker's own
  "position in repeat" columns.
- `reciprocal_overlap` — `min()` of the two above; both sequences must be
  substantially explained for this to be high, so it's the number to sort/
  filter on if you want to distinguish "your motif basically *is* this
  known repeat" from "your motif happens to embed a small fragment of it"
  or "this known repeat happens to be one small piece of your much larger
  motif" (either of those shows up as one of the two individual percentages
  being high while the other is low).

If a query was doubled, a real match that happens to span the artificial
doubling junction can produce an inflated or split-looking raw `sw_score`
right at that boundary — this is an accepted approximation of the doubling
trick that only `pct_query_covered` corrects for (by capping), not
`sw_score`/`pct_divergence`/etc., which are reported exactly as
RepeatMasker output them. Treat a borderline or surprising top hit as worth
a manual look rather than taking the score at face value.

**`ERROR:__main__:FamDB data directory not found`** (from
`run_repeatmasker_known_screen.log`) **is expected on a freshly-built env,
not a bug** — a one-time manual setup step is required every time this
conda env is rebuilt. `workflow/envs/repeatmasker.yaml` is floor-pinned to
the latest release (`repeatmasker>=4.2.4`, same convention as
`mafft>=7.487` elsewhere in this repo), and on that version chain the
Dfam/FamDB library comes from a separate `famdb` conda package that —
confirmed by reading `recipes/famdb/build.sh` in bioconda-recipes — ships
only the `famdb.py` tool, no `*.h5` data at all. `-species` mode never
works out of the box here, on any node, regardless of internet access.

An older pin (`repeatmasker=4.1.5`, matching `vollgerlab/Rhodonite`'s
already-working env on this cluster) was tried instead, since its
`post-link.sh` downloads a library automatically. That was rejected after
checking it directly: it's a frozen, pre-2023-schema single-file snapshot
(curated-only, Dfam release 3.7) with two problems layered on top of each
other — no uncurated/RepeatModeler-derived content (where species-specific
families most likely live), and its own bundled `famdb.py` (v0.4.2)
predates FamDB's schema changes in both v1.0 (Nov 2023) and v3.0.0
(May 2026), so it **cannot read current-format files from dfam.org at
all** even if you wanted to manually upgrade its library. Concretely, for
zebra finch, `famdb.py -i <path-to-Dfam.h5> lineage -a "zebra finch"`
against that 4.1.5-bundled library showed **every node from Aves down to
the species itself at `[0]`** — Aves, Neognathae, Passeriformes,
Passeroidea, Estrildidae, Estrildinae, Taeniopygia, Taeniopygia guttata,
all zero; everything `-species` pulled in came from broad ancestor clades
(Amniota, Vertebrata, Sauropsida) almost certainly dominated by unrelated
model organisms. A clean `known_repeat_hits.tsv` result would have meant
nothing against that library. Floor-pinning to latest costs a manual setup
step but is the only path to current, full (curated+uncurated) data.

**One-time setup** (per conda env build — repeat this if the env is ever
rebuilt: env yaml change, `--conda-cleanup-envs`, a fresh clone, or a
different `--conda-prefix` all invalidate it):
1. Locate this rule's `famdb.py` and confirm its version first —
   `famdb.py`'s own CLI has changed twice (schema v1.0 and v3.0.0); don't
   assume flags from any documentation, including this README, match your
   installed copy without checking `--help` yourself:
   ```
   <conda-env>/share/RepeatMasker/famdb.py info
   ```
2. Current Dfam (as of this pipeline's `repeatmasker>=4.2.4` pin) uses
   FamDB format v3: 4 independently-partitioned components — Curated
   Consensus (`cc`), Curated HMMs (`ch`), Uncurated Consensus (`uc`),
   Uncurated HMMs (`uh`) — plus an always-required root file. From a
   machine with internet access, download at least the root file from
   `https://www.dfam.org/releases/current/families/FamDB/` into an empty
   directory.
3. Ask `famdb.py` itself which partition files you actually need, rather
   than guessing filenames:
   ```
   famdb.py -i <dir> check "<species>"
   ```
   Download the `cc` + `uc` partitions it reports (consensus sequences,
   what RepeatMasker's default rmblast search engine uses; `ch`/`uh` HMM
   variants are only needed for higher-sensitivity nhmmer-based searches).
4. Transfer everything to the cluster if it isn't already there, then
   inside the conda env run that RepeatMasker install's own `configure`
   script pointing at the directory: `-famdb_dir /path/to/downloaded/files`.
   `run_repeatmasker_known_screen` already passes `-uncurated` to
   RepeatMasker so that component actually gets searched, not just curated.

**Before trusting a clean screen as conclusive**, re-run the same lineage
check against whatever you just configured —
`famdb.py -i <dir> lineage -a "<species>"` — and confirm it isn't all
zeros the way the 4.1.5 library was. A `has_known_hit=False` result is only
informative if the library actually had something to compare against.

## `results/summary_table.tsv` columns

| column | meaning |
|---|---|
| `label` | unique ID, `{motif}_rank{rank}_n{N}` — matches the FASTA header in `repeatmasker_custom_lib.fasta` |
| `motif` / `target_period` / `window` | which scan bin this cluster came from |
| `rank` | this cluster's size rank *within its bin* (1 = most loci in that bin) |
| `n_input_sequences` | loci that joined this specific cluster — the core "how much support" number |
| `pct_of_bin` | `n_input_sequences` as a % of the bin's total loci |
| `consensus_length` / `gc_content` | of this cluster's consensus sequence |
| `is_redundant` | True if this cluster lost a redundancy comparison to another (higher-support) cluster |
| `redundant_with` | the winning cluster's `label`, if `is_redundant` |
| `group_size` | how many clusters (across the whole run) got grouped together as redundant with each other |

## Sanity checks worth doing on real data

- Sort `summary_table.tsv` by `n_input_sequences` descending — your
  strongest centromeric/large-satellite candidates are at the top,
  regardless of which bin they came from.
- Check `group_size` for anything > 1 — that cluster was found to be
  redundant with something else; `redundant_with` tells you which cluster
  won and why (compare `n_input_sequences` between the two).
- A bin producing zero clusters (missing from `all_clusters.tsv` entirely)
  means every locus in that bin failed to join a cluster of `min_cluster_size`
  — the bin wasn't coherent enough to be trustworthy, not a bug.
- Eyeball a cluster's `ranked_clusters/rank0N_*.fasta` against
  `01_raw_consensus.fasta` in an alignment viewer before fully trusting a
  borderline consensus, especially one near the `min_cluster_size` floor.
