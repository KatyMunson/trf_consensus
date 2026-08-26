# TRF Motif Consensus-of-Consensus Pipeline

Fully automatic, no manually-curated candidate list: given a **manifest**
of one or more entries (one individual analyzed by one assembly method
each), it scans each entry's raw TRF `.dat` file(s) for candidate period
bins, builds a ranked set of consensus sequences per bin via iterative
clustering, pools every entry's surviving clusters and resolves redundancy
across them (including tandem-duplicate relationships and cross-entry
matches), and outputs a single RepeatMasker-ready custom library plus a QC
summary table annotated with which individuals/methods found each family.
A single-assembly run is just a 1-row manifest — there's no separate
single-`.dat` mode.

## Pipeline

```
config/manifest.tsv
        |
        v  (once per entry, namespaced under results/{entry_id}/...)
scan_dat_candidates (checkpoint)
        |
        v
extract_by_period  (per bin — pools both haplotypes' .dat files for a phased entry)
        |
        v
rank_family_clusters  (per bin: iterative clustering + per-cluster consensus + total_copy_number)
        |
        v
build_all_clusters_table  +  combine_cluster_fastas
        |
        v  (once every entry has finished the above)
plot_copy_number_diagnostic  (results/copy_number_diagnostic.png)
plot_copy_number_qc_diagnostic  (results/copy_number_qc_scatter.png)
        |
        v
filter_and_pool_clusters  (drop below min_total_copy_number, pool across entries)
        |
        v
cross_cluster_comparison  (sharded)  ->  concatenate_cross_cluster_comparison
        |
        v
resolve_redundancy  (names families, annotates provenance)
        |
        v
results/repeatmasker_custom_lib.fasta + results/summary_table.tsv
        |
        v
plot_top_families  (per-individual + global ranking plots)
plot_copy_number_vs_recurrence  (results/copy_number_vs_recurrence.png)
```

1. **`scan_dat_candidates`** scans one entry's whole `.dat` file(s) directly
   (both haplotypes pooled together for a `phased_pooled` entry) and
   nominates every period bin meeting the `candidate_scan` filters (see
   below). This is a Snakemake **checkpoint**, run once per manifest entry —
   that entry's bin list isn't known until this actually runs, so everything
   downstream for that entry is generated dynamically from its output. There's
   no manual pause or curation step: every passing bin proceeds automatically.
2. **`extract_by_period`** pulls each bin's loci out of the entry's `.dat`
   file(s). For a `phased_pooled` entry this is the haplotype-pooling step:
   every matching locus from both haplotype files is written into one
   combined `01_raw_consensus.fasta`/`.tsv`, so a phased and an unphased
   entry both end up representing the same amount of underlying diploid
   genome (no ~2x copy-number gap between them purely from file count).
3. **`rank_family_clusters`** replaces the old two-round primary/secondary
   system entirely. For each bin, it repeatedly: picks a reference sequence
   (median length among whatever remains), anchor-matches everyone else
   against it (rotation- and strand-aware, via `edlib`), peels off what
   matches as one cluster, and repeats on the leftovers — until the
   remaining pool is smaller than `min_cluster_size` or `max_rounds` is
   hit. All discovered clusters in a bin are then sorted by size: rank 1 =
   most support, rank 2 = next, etc. Each cluster gets its own MAFFT
   alignment and majority-rule consensus, plus a summed `total_copy_number`
   (raw TRF `copy_number`, not locus count) across every locus that joined
   it — the per-entry support metric everything downstream ranks and filters
   on — and copy-number-weighted `mean_percent_match` / `mean_entropy`
   (TRF's own per-locus quality fields) across those same loci, feeding
   step 6's QC cross-reference diagnostic. **Sequences that never join a
   cluster of the minimum size are
   dropped, not reported** — this pipeline is tuned for finding high-signal
   candidates (centromeric/large satellite arrays), not for characterizing
   background noise.

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
   every bin's ranked clusters into one master table and one combined FASTA
   — per entry.
5. **`plot_copy_number_diagnostic`** runs once every entry has finished the
   above, pooling every entry's *pre-filter* `total_copy_number` values
   (including clusters about to be dropped — that's the point) into one
   log-x-scaled histogram with a vertical line at the configured
   `min_total_copy_number`, so you can eyeball whether the threshold
   actually sits in a gap between background noise and real signal before
   committing to a rerun.
6. **`plot_copy_number_qc_diagnostic`** runs alongside step 5 (same
   pre-filter input timing), because the plain histogram alone can't
   always resolve an ambiguous noise/signal boundary — a real run might
   show a noisy plateau with no clean valley rather than the two-population
   split the histogram assumes. Cross-references the same
   `total_copy_number` against `target_period` and two of TRF's own
   per-locus quality fields aggregated per cluster
   (`mean_percent_match`, `mean_entropy`, both **copy-number-weighted**
   across the cluster's raw loci — a locus that itself represents a large
   array counts proportionally more toward the cluster's average quality
   than a low-copy-number one), colored by whether each cluster currently
   passes the threshold. Look for red (filtered) points sitting high on
   percent-match/entropy — candidates possibly excluded for a reason other
   than being noise (e.g. an array fragmented across assembly gaps) — and
   blue (passing) points sitting low on those axes, which look questionable
   on every other available signal.
7. **`filter_and_pool_clusters`** drops each entry's clusters below
   `min_total_copy_number`, then pools every surviving cluster from every
   entry into one table and one FASTA. Clusters are renamed
   `cluster_uid = {entry_id}__{source_cluster_label}` so labels that
   collide across entries (e.g. two entries both producing a
   `cand171_rank1_n...` label) stay distinct once pooled.
8. **`cross_cluster_comparison`** compares every pooled cluster found — all
   entries, all bins, all ranks, including different ranks within the same
   bin — using
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

   With multi-entry pooling, the number of pooled clusters compared here
   (and thus the O(N²) pair count) grows well beyond what the earlier
   single-assembly pipeline was sized for, so this step gets three
   performance changes on top of the unchanged comparison logic above:
   - An **exact, lossless pre-filter**: a pair whose length ratio is too
     large for `max_copies` tiled copies to ever reach `min_coverage` is
     skipped before any alignment work, with no accuracy tradeoff — it
     falls straight out of the existing `max_copies`/`min_coverage` values.
   - An **opt-in approximation**, `null_test_skip_margin` (default `0.0`,
     off), that auto-confirms a candidate pair without running the shuffle
     test when it clears both `min_coverage` and `min_identity` by at least
     the margin — see the config section below before enabling it.
   - **Sharding**: the rule runs as `num_shards` parallel Snakemake jobs,
     each processing a contiguous slice of the pair list via
     `--shard-index`/`--num-shards`, then `concatenate_cross_cluster_comparison`
     joins them back into the same `results/cross_cluster_comparison.tsv`
     filename and schema `resolve_redundancy` already expects. Each
     candidate pair's shuffle-null RNG is seeded deterministically from
     `(shuffle_seed, label_A, label_B)` rather than one shared sequential
     stream, so results are identical regardless of shard count or
     execution order.
9. **`resolve_redundancy`** processes pooled clusters in descending
   `total_copy_number` order (not `n_input_sequences` — copy number is the
   fair cross-entry/cross-method support metric now that pooling is
   involved). Each cluster is checked only against clusters **already
   confirmed as winners** — never against another loser, never chained
   through an intermediate. If it's flagged against an existing winner, it
   joins that winner's family; otherwise it becomes a new winner and the
   representative of a new family. **This deliberately avoids transitive
   union-find** (A~B flagged + B~C flagged does NOT imply A~C) — plain
   union-find over flagged pairs was tried first and, on real data,
   collapsed dozens of genuinely unrelated candidates spanning an 11x
   period range into one meaningless "family" simply because each
   consecutive pair along a chain happened to be individually flagged.
   Every winner is named `SAT{motif_length}_{letter}`, where `motif_length`
   is that winner's own exact `consensus_length` (no rounding/averaging —
   different entries' independent MSAs can legitimately disagree by a base
   or two) and `letter` disambiguates families that happen to share a
   motif_length, assigned in discovery order. Every row is annotated with
   which entries/individuals/methods its family was found in.

   **Same-entry consolidation pre-pass**: before that global pass runs, the
   exact same non-chained algorithm runs once per entry, scoped to that
   entry's own clusters only, resolving them down to one "sub-winner" per
   real family present in that entry — only sub-winners are fed into the
   global pass. This exists because a single entry's TRF period scan can
   nominate several candidate bins for the same real monomer (period jitter
   surviving `nms_radius` suppression, e.g. periods 411/417/423 all mutually
   flagged similar) — without this pre-pass, two such same-entry clusters
   could each independently get flagged against the same external winner
   and land in the same family without ever being compared to each other,
   colliding on the RepeatMasker library's `{final_name}__{entry_id}`
   header identity. A cluster consolidated away within its entry records
   which sub-winner it was folded into via `within_entry_consolidated_into`
   (a separate, earlier provenance step from the global `is_representative`
   demotion) and inherits that sub-winner's eventual family. Ranking always
   uses a sub-winner's own `total_copy_number` as-is, never summed across
   what it consolidated — bin extraction windows (`period ± window`) can
   overlap, so summing risks double-counting the same raw TRF loci.

   Nothing is dropped from the summary table — every original cluster call
   still gets a row, non-representative members kept and marked
   `is_representative=False`. The final RepeatMasker library is narrower:
   it gets exactly one sequence per entry per family (that entry's own
   sub-winner, win or lose globally), demoted to
   `#Satellite/redundant_with_<final_name>` in the library's classification
   field when it lost the global pass — a within-entry-consolidated
   cluster doesn't get its own library entry, since it's redundant with its
   own entry's sub-winner sequence by construction.
10. **`plot_top_families`** reads the final summary table and produces one
    ranking bar chart per individual (that individual's own families, by the
    *max* `source_copy_number` across that individual's own entries — not a
    sum, to avoid conflating two methods' independent measurements of the
    same underlying quantity) plus one global chart across all
    individuals/entries. Bars are colored by whether the family was found in
    more than one individual, as a visual cross-check against
    `n_individuals_present`.
11. **`plot_copy_number_vs_recurrence`** runs after step 9, one point per
    family (the representative's own `total_copy_number`) against how many
    entries and how many methods independently confirmed it
    (`n_entries_found` / `n_methods_confirming`, as two side-by-side panels
    since they can tell different stories — e.g. a family found in several
    entries but only one method is a weaker claim than one confirmed by two
    different methods). Cross-entry recurrence is confirmation a
    single-entry copy-number threshold can't see at all; this view is what
    makes a *more permissive* `min_total_copy_number` combined with a good
    recurrence signal a defensible choice, rather than asking one static
    threshold to do all the noise/signal separation alone.

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

## Manifest

`config/manifest.tsv` (see the worked example at
`config/manifest.example.tsv`) is a tab-separated file, one row per
**entry** = one individual analyzed by one assembly method:

| column | meaning |
|---|---|
| `entry_id` | `{individual}_{assembly_method}` — the unit of analysis from here on, and the namespace prefix for that entry's `results/{entry_id}/...` outputs |
| `individual` | which biological individual |
| `assembly_method` | which assembly pipeline produced the input(s) |
| `phasing_status` | `phased_pooled` (two haplotype-specific `.dat` files, pooled at calling time) or `unphased` (one diploid `.dat` file) |
| `fasta_paths` | one path, or two `;`-separated paths (`hap1;hap2`) — **provenance/documentation only**, no script actually reads these |
| `trf_dat_paths` | one path, or two `;`-separated paths, positionally paired with `fasta_paths` — this is what the pipeline actually consumes |

A `phased_pooled` entry must give exactly 2 paths in both `fasta_paths` and
`trf_dat_paths`; an `unphased` entry must give exactly 1 — the Snakefile
validates this at load time and fails fast with a clear error otherwise.
Pooled `.dat` files for one entry should use distinct/prefixed contig names
across the two haplotypes if possible; if they don't,
`extract_by_period.py` defensively prefixes locus ids with a haplotype-file
index (`h0:`, `h1:`) to avoid silent id collisions in per-locus lookups
downstream.

A single-assembly run is just a 1-row manifest.

## Setup

Edit `config/config.yaml`:
- `manifest`: path to your manifest (see "Manifest" above)
- `trf_dat_default_seqname`: only needed if a `.dat` has no header lines at
  all — one global fallback applied to every entry, not configurable per entry
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

**Both `--use-envmodules` and `--use-conda` are required together** — most
rules use envmodules (matching the rest of this pipeline's cluster-module
setup), but a few deliberately have no `envmodules:` fallback:
`run_repeatmasker_known_screen` (it installs its own fresh, pipeline-owned
RepeatMasker via conda rather than reusing whatever's already on the
cluster — see "Known-repeat screening" below), and
`plot_copy_number_diagnostic`/`plot_copy_number_qc_diagnostic`/
`plot_top_families`/`plot_copy_number_vs_recurrence` (matplotlib isn't
assumed to have a cluster module, so they're conda-only). Passing both
flags lets Snakemake use envmodules wherever a rule defines them and fall
back to conda for the rest; passing `--use-envmodules` alone would run
those rules with no environment management at all.

(Swap `--use-envmodules` for `--use-conda` to run off the per-rule envs in
`workflow/envs/` instead: `python_base.yaml` for pure-Python rules,
`edlib.yaml` for `cross_cluster_comparison`, `cluster_rank.yaml` [edlib + mafft together] for `rank_family_clusters`,
`plotting.yaml` [matplotlib] for the four plotting rules
(`plot_copy_number_diagnostic`, `plot_copy_number_qc_diagnostic`,
`plot_top_families`, `plot_copy_number_vs_recurrence`).)

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
permissive (less strict) confirmation criterion. Dropping `n_shuffles` to
roughly 8-10 is a reasonable way to roughly halve the null-test stage's
cost while `max_null_passes` stays at 0, if runtime is still an issue after
sharding — see the comment above the value in `configexample.yaml`.

`num_shards` (default 4) splits the pairwise comparison across this many
parallel Snakemake jobs (see the `cross_cluster_comparison` step above) —
raise it for a larger pooled-cluster count / more available SGE slots, or
set it to 1 to run as a single job like before sharding was added.
`null_test_skip_margin` (default 0.0, off) is the opt-in shuffle-test
shortcut described above — leave it at 0.0 unless you've validated a
nonzero margin against a `margin=0` baseline run on your own data (compare
`n_flagged` in the log between the two), since it changes the statistical
guarantee for confirmed pairs. `configexample.yaml` documents two starting
points (0.05 moderate, 0.10 conservative) if you decide to enable it.

**`copy_number_filter`** — `min_total_copy_number` (default 500) is the
per-entry cluster filter applied in `filter_and_pool_clusters`, before
cross-entry pooling. Check `results/copy_number_diagnostic.png` (a
log-scaled histogram of every entry's pre-filter `total_copy_number`
values with this threshold marked) to confirm it sits in a real gap
between background noise and true arrays before trusting a run — and if
the histogram alone doesn't show a clean gap (a real, not hypothetical,
failure mode: a noisy plateau with no clear valley), cross-check against
`results/copy_number_qc_scatter.png` (copy number vs. period length /
percent-match / entropy) and, after a run completes,
`results/copy_number_vs_recurrence.png` (copy number vs. how many
entries/methods independently confirmed each family) before deciding the
threshold is placed well. A middling-copy-number family confirmed across
several entries/methods is stronger evidence than an equally-scored family
found only once — recurrence a single static per-entry threshold can't
see on its own, and a reason a more permissive threshold can be a
defensible choice once you're leaning on that cross-check too.

**`visualization`** — `top_n_families` (default 20) controls how many bars
`plot_top_families` draws per individual and in the global ranking plot.

**`known_repeat_screen`** — `species` is the RepeatMasker `-species` value;
`double_sequences` (default true) controls the rotation-tolerant doubling
trick. See "Known-repeat screening" below for the full rationale.

## Outputs

Per-entry (namespaced under `results/{entry_id}/...`):
- `results/{entry_id}/candidate_periods.tsv` / `candidate_periods_raw.tsv` /
  `candidate_periods_manifest.tsv` — that entry's scan output (see
  `scan_dat_candidates.py` docstring for column definitions)
- `results/{entry_id}/{motif}/01_raw_consensus.fasta` + `.tsv` — extracted
  per-locus consensus motifs + provenance for that bin (haplotype-pooled
  for a `phased_pooled` entry)
- `results/{entry_id}/{motif}/ranked_clusters/rank01_n<N>.fasta`,
  `rank02_...`, etc. — one consensus per cluster found in that bin
- `results/{entry_id}/{motif}/ranked_clusters_summary.tsv` — per-bin
  cluster ranking, including `total_copy_number` and the
  copy-number-weighted `mean_percent_match` / `mean_entropy` (TRF's own
  per-locus quality fields, averaged across the cluster's raw loci)
- `results/{entry_id}/all_clusters.tsv` — every cluster from every bin of
  this entry, one master table
- `results/{entry_id}/all_clusters_consensus.fasta` — every cluster's
  consensus for this entry, combined

Across entries:
- `results/copy_number_diagnostic.png` — pre-filter `total_copy_number`
  distribution across every entry, with the `min_total_copy_number`
  threshold marked
- `results/copy_number_qc_scatter.png` — the same pre-filter clusters,
  `total_copy_number` (log x-axis, colored by pass/fail) vs. `target_period`
  / `mean_percent_match` / `mean_entropy`, for when the histogram alone
  doesn't show a clean noise/signal gap
- `results/filtered_pooled_clusters.tsv` — every entry's surviving clusters
  (post-filter), pooled, `cluster_uid = {entry_id}__{source_cluster_label}`
- `results/pooled_consensus.fasta` — every surviving cluster's consensus,
  headers = `cluster_uid`
- `results/cross_cluster_comparison.tsv` — pairwise similarity/multiple-of
  check across every pooled cluster
- `results/repeatmasker_custom_lib.fasta` — **final output**: every
  pooled cluster's consensus, headers as
  `>{final_name}__{entry_id}#{classification}` (family representatives) or
  `>{final_name}__{entry_id}#{classification}/redundant_with_{final_name}`
  (other family members)
- `results/summary_table.tsv` — **final output**: one row per original
  cluster call, named/grouped into families with provenance — see the
  columns table below
- `results/known_repeat_hits.tsv` — **final output**: every library entry
  joined against a Dfam/RepBase known-repeat screen, see "Known-repeat
  screening" below
- `results/top_families_{individual}.png`, `results/top_families_global.png`
  — **final output**: family-ranking bar charts, see "Pipeline" step 10
- `results/copy_number_vs_recurrence.png` — **final output**: one point per
  family, representative `total_copy_number` vs. `n_entries_found` /
  `n_methods_confirming`, see "Pipeline" step 11

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

`run_repeatmasker_known_screen` runs RepeatMasker from inside
`results/known_repeat_screen/`, not the repo root — RepeatMasker's own
`RM_<pid>.<timestamp>` scratch directory is always created relative to the
process's actual working directory (confirmed by reading its source,
`createTempDir()`; `-dir` is only ever consulted as a fallback if writing
to cwd fails outright), so without this it would litter the repo root on
every invocation. RepeatMasker does normally self-delete its own scratch
dir, but only via the very last line of its main script — any crash or
interruption anywhere earlier skips that entirely, which is why you may
still see a stray `RM_*` directory after a failed run despite this. The
rule also traps on exit to remove it regardless of success or failure, so
this should be rare going forward; any that are still lying around from
before this fix are safe to delete by hand
(`rm -rf results/known_repeat_screen/RM_*` — or, from old runs before this
was contained to that directory, `rm -rf RM_*` in the repo root).

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

**Reading `known_repeat_hits.tsv`**: every library entry (identified by its
FASTA header up to the first `#`, i.e. `{final_name}__{entry_id}` — the
same join key you'd reconstruct from `summary_table.tsv`'s `final_name` +
`source_entry_id`) is present, but **this is not guaranteed to be one row
per entry** — an entry with no hit gets exactly one row
(`has_known_hit=False`, rest `NA`), a normal, common, and often *expected*
outcome, not a failure (much of the satellite DNA in a non-model species
genome is genuinely undescribed in existing databases); an entry with hits
gets **one row per distinct matched repeat name**, so a motif that matches
two unrelated known families shows up as two rows. If you need exactly one
row per entry (e.g. for a simple join), sort by `reciprocal_overlap`
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

One row per original cluster call — nothing dropped, non-representative
family members are kept and marked `is_representative=False`. Note that
`is_representative=False` covers two distinct demotions now (see "Pipeline"
step 9's same-entry consolidation pre-pass): a row can be demoted globally
(lost the cross-entry pass to some other entry's stronger cluster) and/or
consolidated within its own entry first (folded into that entry's own
sub-winner before the global pass ever ran) — `within_entry_consolidated_into`
is what distinguishes the latter.

| column | meaning |
|---|---|
| `final_name` | `SAT{motif_length}_{letter}` — the harmonized family name. Matches the FASTA header prefix in `repeatmasker_custom_lib.fasta` (`{final_name}__{source_entry_id}`) |
| `motif_length` | the family's representative cluster's own exact `consensus_length` (no rounding/averaging) |
| `is_representative` | True for the single row per family with the highest `source_copy_number` overall — its consensus is what goes in the RepeatMasker library unqualified; every other row (globally demoted or within-entry-consolidated) is `False` |
| `within_entry_consolidated_into` | `NA` if this row was its own entry's sub-winner (went into the global pass directly); otherwise the `source_cluster_label` of the sub-winner it was folded into within its own entry, before the global pass ran |
| `source_entry_id` / `source_individual` / `source_assembly_method` / `source_phasing_status` | which manifest entry this specific cluster call came from |
| `source_cluster_label` | that entry's own cluster label (`{motif}_rank{rank}_n{N}`), before pooling |
| `source_copy_number` | this cluster's own `total_copy_number` (summed raw TRF copy number, not locus count) — the ranking key resolve_redundancy used at both the within-entry and global level (never summed across consolidated duplicates) |
| `source_consensus_length` / `source_gc_content` | of this specific cluster's own consensus (can differ slightly from `motif_length` — see "Pipeline" step 9) |
| `n_entries_found` / `entries_found` | how many / which manifest entries this family was found in at all |
| `n_individuals_present` / `individuals_found` | how many / which individuals this family was found in |
| `n_methods_confirming` / `methods_found` | how many / which assembly methods this family was found in |

## Sanity checks worth doing on real data

- Sort by `is_representative` then `source_copy_number` descending within a
  family, or just eyeball `results/top_families_global.png` — your
  strongest, most-confirmed candidates are the ones with a tall bar colored
  "found in >1 individual."
- `n_individuals_present` / `n_methods_confirming` are the new cross-checks
  worth scanning: a family private to one individual *and* one method is a
  weaker claim than one confirmed across individuals and methods — compare
  against `results/copy_number_diagnostic.png` too, since a private,
  low-copy-number family sitting just above the filter threshold is worth
  a second look.
- `SAT{motif_length}_a` / `SAT{motif_length}_b` sharing a `motif_length` are
  **not** the same family — the letter, not the length, disambiguates them;
  don't assume they should be merged just because they share a number.
- A bin producing zero clusters (missing from that entry's `all_clusters.tsv`
  entirely) means every locus in that bin failed to join a cluster of
  `min_cluster_size` — the bin wasn't coherent enough to be trustworthy, not
  a bug.
- Eyeball a cluster's `ranked_clusters/rank0N_*.fasta` against
  `01_raw_consensus.fasta` in an alignment viewer before fully trusting a
  borderline consensus, especially one near the `min_cluster_size` floor.
- Ties in `source_copy_number` when resolving which cluster becomes a
  family's representative break deterministically but not obviously — via
  the order clusters appear in `filtered_pooled_clusters.tsv` (itself
  manifest-row order × per-entry bin/rank order), not an explicit tiebreak
  rule. Not usually worth worrying about, but worth knowing if two runs
  ever appear to disagree on which near-identical-support cluster won.
