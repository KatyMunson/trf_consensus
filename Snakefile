# =============================================================================
# Snakefile — TRF Motif Consensus-of-Consensus Pipeline (multi-entry, fully
# programmatic)
#
# Manifest-driven: every run is defined by config["manifest"] (see
# config/manifest.example.tsv), one row per "entry" = one individual
# analyzed by one assembly method. A phased entry pools both haplotypes'
# TRF .dat files at calling time (before any clustering); an unphased entry
# has one .dat file. A single-assembly run is just a 1-row manifest — there
# is no separate single-.dat mode.
#
# Per entry (namespaced under results/{entry_id}/...), unchanged from the
# single-assembly pipeline in its core logic:
#
#   1. scan_dat_candidates (checkpoint) — scans the entry's .dat file(s)
#      directly and nominates every period bin meeting the candidate_scan
#      filters. No manual pause/curation step — every passing bin proceeds.
#   2. extract_by_period — per bin, pull consensus motifs (TRF .dat field 14)
#      for entries within period +/- window, pooling every matching locus
#      across all of this entry's .dat file(s).
#   3. rank_family_clusters — per bin, iteratively peel off coherent
#      sequence clusters (rotation/strand-aware anchor matching) and rank
#      them by size, also summing raw TRF copy_number per cluster into
#      total_copy_number (the per-entry support metric used downstream).
#   4. build_all_clusters_table + combine_cluster_fastas — aggregate every
#      bin's ranked clusters into one per-entry master table/FASTA.
#
# Across entries:
#
#   2b. plot_copy_number_diagnostic — once every entry's Stage 2 has run,
#       pools every entry's PRE-FILTER total_copy_number values into one
#       log-scaled histogram with the configured min_total_copy_number
#       threshold marked, so it can be sanity-checked before committing to
#       a rerun.
#   3.  filter_and_pool_clusters — drops each entry's clusters below
#       min_total_copy_number, then pools every survivor across every entry
#       into one table (cluster_uid = {entry_id}__{source_cluster_label})
#       and one concatenated FASTA.
#   4.  cross_cluster_comparison — pairwise rotation/strand-aware comparison
#       across every pooled cluster (all entries, all bins, all ranks),
#       flagging near-identical or integer-multiple pairs. Unchanged logic —
#       operates purely on sequence content.
#   5.  resolve_redundancy — processes pooled clusters in descending
#       total_copy_number order; each cluster is checked only against
#       already-confirmed winners (never against another loser, never
#       chained transitively through an intermediate). Winners become named
#       families (SAT{motif_length}_{letter}); other members are marked
#       redundant with their family, not dropped. Every row is annotated
#       with which entries/individuals/methods it was found in.
#   6.  plot_top_families — per-individual and global bar charts ranking
#       families by copy number (max across an individual's own entries,
#       per spec — not summed, to avoid conflating two methods'
#       measurements of the same underlying quantity).
#   7.  known-repeat screening (prepare_known_repeat_query ->
#       run_repeatmasker_known_screen -> parse_known_repeat_hits) — runs
#       RepeatMasker -species against OUR OWN discovered consensus motifs
#       (not the assembly) to check whether any match a previously
#       characterized repeat family in Dfam/RepBase. Unchanged — already
#       treats FASTA headers generically. Uses a fresh, pipeline-owned
#       RepeatMasker install (conda only, see workflow/envs/repeatmasker.yaml
#       — floor-pinned to the latest release. The FamDB library it needs
#       does NOT come bundled or auto-downloaded on this version chain, so
#       the first run after a fresh env build will fail with "FamDB data
#       directory not found" — this is expected, not a bug; see the
#       README's "Known-repeat screening" section for the one-time manual
#       setup), independent of whatever version/library may already be on
#       the cluster.
#
# This requires Snakemake checkpoints (checkpoint scan_dat_candidates): each
# entry's bin list isn't known until its own scan actually runs, so
# downstream rules use input functions that call
# checkpoints.scan_dat_candidates.get(entry_id=...) to discover it at
# DAG-re-evaluation time, rather than a static manifest.
#
# Usage: snakemake -s Snakefile --configfile config/config.yaml --cores <N> \
#          --use-envmodules --retries 3
# See README.md for full documentation and cluster submission instructions.
#
# Author:  KM
# Created: 2026-08
# =============================================================================

import os

wildcard_constraints:
    entry_id="[^/]+",
    motif="[^/]+",


# --- manifest parsing (static, at load time -- NOT a checkpoint) ---
# First non-blank line is a mandatory header row (entry_id, individual,
# assembly_method, phasing_status, fasta_paths, trf_dat_paths) and is
# skipped, matching config/manifest.example.tsv's format.
MANIFEST = {}
with open(config["manifest"]) as f:
    header_skipped = False
    for lineno, line in enumerate(f, 1):
        if not line.strip():
            continue
        if not header_skipped:
            header_skipped = True
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) != 6:
            raise ValueError(
                f"manifest line {lineno}: expected 6 tab-separated fields "
                f"(entry_id, individual, assembly_method, phasing_status, "
                f"fasta_paths, trf_dat_paths), got {len(fields)}"
            )
        entry_id, individual, assembly_method, phasing_status, fasta_paths, trf_dat_paths = fields
        if entry_id in MANIFEST:
            raise ValueError(f"manifest line {lineno}: duplicate entry_id '{entry_id}'")
        if phasing_status not in ("phased_pooled", "unphased"):
            raise ValueError(
                f"manifest line {lineno}: phasing_status must be 'phased_pooled' or "
                f"'unphased', got '{phasing_status}'"
            )
        n_fasta = len(fasta_paths.split(";"))
        n_dat = len(trf_dat_paths.split(";"))
        expected_n = 2 if phasing_status == "phased_pooled" else 1
        if n_fasta != expected_n or n_dat != expected_n:
            raise ValueError(
                f"manifest line {lineno}: phasing_status '{phasing_status}' requires "
                f"exactly {expected_n} ;-separated path(s) in both fasta_paths and "
                f"trf_dat_paths, got {n_fasta} fasta path(s) and {n_dat} trf_dat path(s)"
            )
        MANIFEST[entry_id] = {
            "individual": individual, "assembly_method": assembly_method,
            "phasing_status": phasing_status, "fasta_paths": fasta_paths,
            "trf_dat_paths": trf_dat_paths,
        }

if not MANIFEST:
    raise ValueError(f"No entries found in manifest {config['manifest']}!")

ENTRY_IDS = list(MANIFEST)
INDIVIDUALS = sorted({v["individual"] for v in MANIFEST.values()})


def entry_dat_paths(wildcards):
    return MANIFEST[wildcards.entry_id]["trf_dat_paths"].split(";")


def get_period_manifest(wildcards):
    """Parse one entry's scan checkpoint manifest into {name: (period, window)}.
    Every rule/function that needs that entry's bin list or a bin's
    period/window calls this, which forces the checkpoint to complete first."""
    path = checkpoints.scan_dat_candidates.get(entry_id=wildcards.entry_id).output.manifest
    entries = {}
    with open(path) as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                fields = line.rstrip("\n").split("\t")
                entries[fields[0]] = (int(fields[1]), int(fields[2]))
    return entries


def all_motifs(wildcards):
    return list(get_period_manifest(wildcards).keys())


def motif_period(wildcards):
    return get_period_manifest(wildcards)[wildcards.motif][0]


def motif_window(wildcards):
    return get_period_manifest(wildcards)[wildcards.motif][1]


def all_motif_periods(wildcards):
    m = get_period_manifest(wildcards)
    return [m[name][0] for name in m]


def all_motif_windows(wildcards):
    m = get_period_manifest(wildcards)
    return [m[name][1] for name in m]


rule all:
    input:
        "results/repeatmasker_custom_lib.fasta",
        "results/summary_table.tsv",
        "results/known_repeat_hits.tsv",
        "results/copy_number_diagnostic.png",
        "results/copy_number_qc_scatter.png",
        "results/copy_number_vs_recurrence.png",
        "results/top_families_global.png",
        expand("results/top_families_{individual}.png", individual=INDIVIDUALS)


checkpoint scan_dat_candidates:
    # Scans this entry's whole .dat file(s) directly (pooled together if
    # there are two, one per haplotype); the bin list downstream rules use
    # for this entry is determined here, at runtime, not known when the
    # DAG is built.
    input:
        dat=entry_dat_paths
    output:
        tsv="results/{entry_id}/candidate_periods.tsv",
        raw_tsv="results/{entry_id}/candidate_periods_raw.tsv",
        manifest="results/{entry_id}/candidate_periods_manifest.tsv"
    log:
        "results/logs/{entry_id}/scan_dat_candidates.log"
    threads: config["resources"]["scan_dat_candidates"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["scan_dat_candidates"]["mem"] * attempt,
        hrs=config["resources"]["scan_dat_candidates"]["hrs"]
    conda:
        "workflow/envs/python_base.yaml"
    envmodules:
        "python/3.11"
    shell:
        "python workflow/scripts/scan_dat_candidates.py "
        "--dat {input.dat} "
        "--min-copy-number {config[candidate_scan][min_copy_number]} "
        "--min-blocks {config[candidate_scan][min_blocks]} "
        "--min-single-block-copy-number {config[candidate_scan][min_single_block_copy_number]} "
        "--min-period-length {config[candidate_scan][min_period_length]} "
        "--window {config[candidate_scan][window]} "
        "--nms-radius {config[candidate_scan][nms_radius]} "
        "--out-tsv {output.tsv} --out-raw-tsv {output.raw_tsv} --out-manifest {output.manifest} "
        "--manifest-name-prefix cand > {log} 2>&1"


rule extract_by_period:
    input:
        dat=entry_dat_paths,
        manifest=lambda wc: checkpoints.scan_dat_candidates.get(entry_id=wc.entry_id).output.manifest
    output:
        fasta="results/{entry_id}/{motif}/01_raw_consensus.fasta",
        tsv="results/{entry_id}/{motif}/01_raw_consensus.tsv"
    params:
        period=motif_period,
        window=motif_window
    log:
        "results/logs/{entry_id}/{motif}/extract_by_period.log"
    threads: config["resources"]["extract_by_period"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["extract_by_period"]["mem"] * attempt,
        hrs=config["resources"]["extract_by_period"]["hrs"]
    conda:
        "workflow/envs/python_base.yaml"
    envmodules:
        "python/3.11"
    shell:
        "python workflow/scripts/extract_by_period.py "
        "--dat {input.dat} --period {params.period} --window {params.window} "
        "--out-fasta {output.fasta} --out-tsv {output.tsv} "
        "--default-seqname '{config[trf_dat_default_seqname]}' > {log} 2>&1"


rule rank_family_clusters:
    # Iteratively peels off coherent clusters (rotation/strand-aware anchor
    # matching) from this bin and ranks them by size. Self-contained: calls
    # mafft and does majority-consensus calling internally, since the
    # cluster count per bin isn't known ahead of time.
    input:
        fasta="results/{entry_id}/{motif}/01_raw_consensus.fasta",
        tsv="results/{entry_id}/{motif}/01_raw_consensus.tsv"
    output:
        clusters_dir=directory("results/{entry_id}/{motif}/ranked_clusters"),
        summary="results/{entry_id}/{motif}/ranked_clusters_summary.tsv"
    params:
        motif_name=lambda wc: wc.motif
    log:
        "results/logs/{entry_id}/{motif}/rank_family_clusters.log"
    threads: config["resources"]["rank_family_clusters"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["rank_family_clusters"]["mem"] * attempt,
        hrs=config["resources"]["rank_family_clusters"]["hrs"]
    conda:
        "workflow/envs/cluster_rank.yaml"
    envmodules:
        "python/3.11", "mafft/7.487"  # verify with `module avail mafft` on liger and adjust
    shell:
        "python workflow/scripts/rank_family_clusters.py "
        "--in-fasta {input.fasta} --in-tsv {input.tsv} --out-dir {output.clusters_dir} "
        "--out-summary {output.summary} --motif-name {params.motif_name} "
        "--min-identity {config[cluster_ranking][min_cluster_identity]} "
        "--max-gap-fraction {config[cluster_ranking][max_gap_fraction]} "
        "--min-cluster-size {config[cluster_ranking][min_cluster_size]} "
        "--max-rounds {config[cluster_ranking][max_rounds]} > {log} 2>&1"


rule build_all_clusters_table:
    input:
        summaries=lambda wc: expand("results/{entry_id}/{motif}/ranked_clusters_summary.tsv",
                                     entry_id=wc.entry_id, motif=all_motifs(wc)),
        dirs=lambda wc: expand("results/{entry_id}/{motif}/ranked_clusters",
                                entry_id=wc.entry_id, motif=all_motifs(wc))
    output:
        "results/{entry_id}/all_clusters.tsv"
    params:
        motifs=lambda wc: all_motifs(wc),
        periods=lambda wc: all_motif_periods(wc),
        windows=lambda wc: all_motif_windows(wc)
    log:
        "results/logs/{entry_id}/build_all_clusters_table.log"
    threads: config["resources"]["build_all_clusters_table"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["build_all_clusters_table"]["mem"] * attempt,
        hrs=config["resources"]["build_all_clusters_table"]["hrs"]
    conda:
        "workflow/envs/python_base.yaml"
    envmodules:
        "python/3.11"
    shell:
        "python workflow/scripts/build_all_clusters_table.py "
        "--motifs {params.motifs} --periods {params.periods} --windows {params.windows} "
        "--summaries {input.summaries} --out-tsv {output} > {log} 2>&1"


rule combine_cluster_fastas:
    # Plain concatenation of every cluster's single-record consensus FASTA
    # across every bin of this entry — safe at runtime since Snakemake
    # guarantees the `directory()` outputs from rank_family_clusters are
    # fully materialized before this rule executes (both are declared inputs).
    input:
        summaries=lambda wc: expand("results/{entry_id}/{motif}/ranked_clusters_summary.tsv",
                                     entry_id=wc.entry_id, motif=all_motifs(wc)),
        dirs=lambda wc: expand("results/{entry_id}/{motif}/ranked_clusters",
                                entry_id=wc.entry_id, motif=all_motifs(wc))
    output:
        "results/{entry_id}/all_clusters_consensus.fasta"
    log:
        "results/logs/{entry_id}/combine_cluster_fastas.log"
    threads: config["resources"]["combine_cluster_fastas"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["combine_cluster_fastas"]["mem"] * attempt,
        hrs=config["resources"]["combine_cluster_fastas"]["hrs"]
    shell:
        # nullglob so an empty ranked_clusters/ (a bin that produced zero
        # clusters — expected, not an error) doesn't try to cat a literal
        # non-matching "*.fasta" and trip `set -e`
        "shopt -s nullglob; "
        "for d in {input.dirs}; do "
        "  for f in \"$d\"/*.fasta; do cat \"$f\"; done; "
        "done > {output} 2> {log}"


rule plot_copy_number_diagnostic:
    # Runs once every entry's Stage 2 is complete, pooling every entry's
    # PRE-FILTER total_copy_number values so the min_total_copy_number
    # threshold used in filter_and_pool_clusters can be sanity-checked
    # before committing to a rerun.
    input:
        clusters_tsv=expand("results/{entry_id}/all_clusters.tsv", entry_id=ENTRY_IDS)
    output:
        "results/copy_number_diagnostic.png"
    log:
        "results/logs/plot_copy_number_diagnostic.log"
    threads: config["resources"]["plot_copy_number_diagnostic"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["plot_copy_number_diagnostic"]["mem"] * attempt,
        hrs=config["resources"]["plot_copy_number_diagnostic"]["hrs"]
    conda:
        "workflow/envs/plotting.yaml"
    shell:
        "python workflow/scripts/plot_copy_number_diagnostic.py "
        "--clusters-tsv {input.clusters_tsv} "
        "--min-total-copy-number {config[copy_number_filter][min_total_copy_number]} "
        "--out-png {output} > {log} 2>&1"


rule plot_copy_number_qc_diagnostic:
    # Same input timing as plot_copy_number_diagnostic (every entry's
    # pre-filter all_clusters.tsv) but cross-references total_copy_number
    # against target_period and TRF's own per-locus quality fields
    # (mean_percent_match, mean_entropy), since the plain histogram alone
    # can't resolve an ambiguous noise/signal middle on its own.
    input:
        clusters_tsv=expand("results/{entry_id}/all_clusters.tsv", entry_id=ENTRY_IDS)
    output:
        "results/copy_number_qc_scatter.png"
    log:
        "results/logs/plot_copy_number_qc_diagnostic.log"
    threads: config["resources"]["plot_copy_number_qc_diagnostic"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["plot_copy_number_qc_diagnostic"]["mem"] * attempt,
        hrs=config["resources"]["plot_copy_number_qc_diagnostic"]["hrs"]
    conda:
        "workflow/envs/plotting.yaml"
    shell:
        "python workflow/scripts/plot_copy_number_qc_diagnostic.py "
        "--clusters-tsv {input.clusters_tsv} "
        "--min-total-copy-number {config[copy_number_filter][min_total_copy_number]} "
        "--out-png {output} > {log} 2>&1"


rule filter_and_pool_clusters:
    input:
        clusters_tsv=expand("results/{entry_id}/all_clusters.tsv", entry_id=ENTRY_IDS),
        fastas=expand("results/{entry_id}/all_clusters_consensus.fasta", entry_id=ENTRY_IDS)
    output:
        tsv="results/filtered_pooled_clusters.tsv",
        fasta="results/pooled_consensus.fasta"
    params:
        entry_ids=ENTRY_IDS,
        individuals=[MANIFEST[e]["individual"] for e in ENTRY_IDS],
        assembly_methods=[MANIFEST[e]["assembly_method"] for e in ENTRY_IDS],
        phasing_statuses=[MANIFEST[e]["phasing_status"] for e in ENTRY_IDS]
    log:
        "results/logs/filter_and_pool_clusters.log"
    threads: config["resources"]["filter_and_pool_clusters"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["filter_and_pool_clusters"]["mem"] * attempt,
        hrs=config["resources"]["filter_and_pool_clusters"]["hrs"]
    conda:
        "workflow/envs/python_base.yaml"
    envmodules:
        "python/3.11"
    shell:
        "python workflow/scripts/filter_and_pool_clusters.py "
        "--clusters-tsv {input.clusters_tsv} --fastas {input.fastas} "
        "--entry-ids {params.entry_ids} --individuals {params.individuals} "
        "--assembly-methods {params.assembly_methods} --phasing-statuses {params.phasing_statuses} "
        "--min-total-copy-number {config[copy_number_filter][min_total_copy_number]} "
        "--out-tsv {output.tsv} --out-fasta {output.fasta} > {log} 2>&1"


rule cross_cluster_comparison:
    input:
        fasta="results/pooled_consensus.fasta"
    output:
        "results/cross_cluster_comparison.tsv"
    log:
        "results/logs/cross_cluster_comparison.log"
    threads: config["resources"]["cross_cluster_comparison"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["cross_cluster_comparison"]["mem"] * attempt,
        hrs=config["resources"]["cross_cluster_comparison"]["hrs"]
    conda:
        "workflow/envs/edlib.yaml"
    envmodules:
        "python/3.11"
    shell:
        "python workflow/scripts/cross_motif_comparison.py "
        "--fastas {input.fasta} --out-tsv {output} "
        "--min-identity {config[cross_cluster_comparison][min_identity]} "
        "--min-coverage {config[cross_cluster_comparison][min_coverage]} "
        "--max-copies {config[cross_cluster_comparison][max_copies]} "
        "--n-shuffles {config[cross_cluster_comparison][n_shuffles]} "
        "--max-null-passes {config[cross_cluster_comparison][max_null_passes]} "
        "> {log} 2>&1"


rule resolve_redundancy:
    input:
        clusters_tsv="results/filtered_pooled_clusters.tsv",
        comparison_tsv="results/cross_cluster_comparison.tsv",
        consensus_fasta="results/pooled_consensus.fasta"
    output:
        lib_fasta="results/repeatmasker_custom_lib.fasta",
        summary="results/summary_table.tsv"
    log:
        "results/logs/resolve_redundancy.log"
    threads: config["resources"]["resolve_redundancy"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["resolve_redundancy"]["mem"] * attempt,
        hrs=config["resources"]["resolve_redundancy"]["hrs"]
    conda:
        "workflow/envs/python_base.yaml"
    envmodules:
        "python/3.11"
    shell:
        "python workflow/scripts/resolve_redundancy.py "
        "--clusters-tsv {input.clusters_tsv} --comparison-tsv {input.comparison_tsv} "
        "--consensus-fasta {input.consensus_fasta} "
        "--classification {config[repeatmasker_classification]} "
        "--out-lib-fasta {output.lib_fasta} --out-summary {output.summary} > {log} 2>&1"


rule plot_top_families:
    input:
        summary="results/summary_table.tsv"
    output:
        global_png="results/top_families_global.png",
        per_individual=expand("results/top_families_{individual}.png", individual=INDIVIDUALS)
    params:
        individuals=INDIVIDUALS,
        top_n=config["visualization"]["top_n_families"]
    log:
        "results/logs/plot_top_families.log"
    threads: config["resources"]["plot_top_families"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["plot_top_families"]["mem"] * attempt,
        hrs=config["resources"]["plot_top_families"]["hrs"]
    conda:
        "workflow/envs/plotting.yaml"
    shell:
        "python workflow/scripts/plot_top_families.py "
        "--summary-tsv {input.summary} --individuals {params.individuals} "
        "--out-per-individual {output.per_individual} --out-global {output.global_png} "
        "--top-n {params.top_n} > {log} 2>&1"


rule plot_copy_number_vs_recurrence:
    # Post-harmonization diagnostic: one point per family (representative's
    # own total_copy_number) vs. how many entries/methods independently
    # confirmed it — recurrence a single-entry copy-number threshold can't
    # see at all. Purely descriptive; no filtering/ranking logic here.
    input:
        summary="results/summary_table.tsv"
    output:
        "results/copy_number_vs_recurrence.png"
    log:
        "results/logs/plot_copy_number_vs_recurrence.log"
    threads: config["resources"]["plot_copy_number_vs_recurrence"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["plot_copy_number_vs_recurrence"]["mem"] * attempt,
        hrs=config["resources"]["plot_copy_number_vs_recurrence"]["hrs"]
    conda:
        "workflow/envs/plotting.yaml"
    shell:
        "python workflow/scripts/plot_copy_number_vs_recurrence.py "
        "--summary-tsv {input.summary} --out-png {output} > {log} 2>&1"


# --- Step 7: known-repeat screening ---
# Screens our own discovered consensus motifs (not the assembly) against a
# Dfam/RepBase library via RepeatMasker -species, to check whether any
# match a previously characterized repeat family. Unchanged from the
# single-assembly pipeline: both scripts treat FASTA headers generically
# (split on the first "#"), with no assumption about label format, so the
# new {final_name}__{entry_id}#{classification} headers need zero changes
# here.

rule prepare_known_repeat_query:
    # Strips our own #classification suffix down to the bare label (our
    # own metadata, not a valid RepeatMasker query name, and "#" has
    # special meaning in RepeatMasker's own library format) and optionally
    # doubles each sequence so a database entry starting at a different
    # rotation phase can still align end-to-end.
    input:
        "results/repeatmasker_custom_lib.fasta"
    output:
        "results/known_repeat_screen/query.fasta"
    params:
        double_flag=lambda wc: "--double" if config["known_repeat_screen"]["double_sequences"] else ""
    log:
        "results/logs/prepare_known_repeat_query.log"
    threads: config["resources"]["prepare_known_repeat_query"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["prepare_known_repeat_query"]["mem"] * attempt,
        hrs=config["resources"]["prepare_known_repeat_query"]["hrs"]
    conda:
        "workflow/envs/python_base.yaml"
    envmodules:
        "python/3.11"
    shell:
        "python workflow/scripts/prepare_known_repeat_query.py "
        "--in-fasta {input} --out-fasta {output} {params.double_flag} > {log} 2>&1"


rule run_repeatmasker_known_screen:
    # Deliberately conda-only, no envmodules fallback — unlike every other
    # rule in this pipeline. The point of this rule is a fresh,
    # pipeline-owned RepeatMasker install, its own bundled library, not
    # reuse of whatever version/library happens to already be on the
    # cluster. Requires a one-time manual FamDB download-and-configure step
    # per env build — see workflow/envs/repeatmasker.yaml and the README's
    # "Known-repeat screening" section.
    input:
        "results/known_repeat_screen/query.fasta"
    output:
        "results/known_repeat_screen/query.fasta.out"
    params:
        species=config["known_repeat_screen"]["species"],
        outdir="results/known_repeat_screen",
        input_abs=lambda wc, input: os.path.abspath(input[0])
    log:
        "results/logs/run_repeatmasker_known_screen.log"
    threads: config["resources"]["run_repeatmasker_known_screen"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["run_repeatmasker_known_screen"]["mem"] * attempt,
        hrs=config["resources"]["run_repeatmasker_known_screen"]["hrs"]
    conda:
        "workflow/envs/repeatmasker.yaml"
    shell:
        # RepeatMasker can decline to write a .out file at all when nothing
        # anywhere in the query matches anything in the library (a
        # completely plausible outcome for undescribed satellite DNA) —
        # RepeatMasker still exits 0 in that case, it just doesn't create
        # the file. `&& touch {output}` guarantees the declared output
        # exists whenever RepeatMasker itself succeeded, so Snakemake
        # doesn't fail the job on a legitimate zero-hit result (which
        # parse_repeatmasker_hits.py treats as zero hits, not an error) —
        # but a genuine RepeatMasker crash (non-zero exit) still fails the
        # job normally, since `touch` is never reached.
        # -uncurated: per RepeatMasker's own help text, it searches
        # curated Dfam families only by default. The whole point of this
        # screening step is maximum sensitivity to any prior
        # characterization — curated or not (uncurated/RepeatModeler-
        # derived families are where species-specific content most likely
        # lives) — so there's no reason to withhold it here the way you
        # might for whole-genome masking.
        # `cd {params.outdir}` + trap: RepeatMasker's own source
        # (createTempDir()) creates its RM_<pid>.<timestamp> scratch
        # directory relative to the process's actual cwd, always — `-dir`
        # is only ever consulted as a fallback if writing to cwd fails
        # outright, so without this it litters the repo root (Snakemake
        # jobs run with cwd = the Snakefile's directory) on every
        # invocation. RepeatMasker does self-delete it, but only via the
        # very last line of its main script (`rm -R $tempdir unless
        # $DEBUG`) — any crash/die anywhere earlier skips that entirely.
        # cd'ing into -dir first at least contains any leftover debris
        # inside results/known_repeat_screen/ instead of the repo root; the
        # trap goes further and guarantees cleanup here even when
        # RepeatMasker itself fails, which its own logic does not.
        "(cd {params.outdir} && trap 'rm -rf RM_*' EXIT && "
        "RepeatMasker -species '{params.species}' -pa {threads} -uncurated "
        "-dir . {params.input_abs}) > {log} 2>&1 && touch {output}"


rule parse_known_repeat_hits:
    input:
        rm_out="results/known_repeat_screen/query.fasta.out",
        query_fasta="results/known_repeat_screen/query.fasta"
    output:
        "results/known_repeat_hits.tsv"
    params:
        doubled_flag=lambda wc: "--doubled" if config["known_repeat_screen"]["double_sequences"] else ""
    log:
        "results/logs/parse_known_repeat_hits.log"
    threads: config["resources"]["parse_known_repeat_hits"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["parse_known_repeat_hits"]["mem"] * attempt,
        hrs=config["resources"]["parse_known_repeat_hits"]["hrs"]
    conda:
        "workflow/envs/python_base.yaml"
    envmodules:
        "python/3.11"
    shell:
        "python workflow/scripts/parse_repeatmasker_hits.py "
        "--rm-out {input.rm_out} --query-fasta {input.query_fasta} --out-tsv {output} "
        "{params.doubled_flag} > {log} 2>&1"
