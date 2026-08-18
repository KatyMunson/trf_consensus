# =============================================================================
# Snakefile — TRF Motif Consensus-of-Consensus Pipeline (fully programmatic)
#
# Fully automatic, no manually-curated periods.tsv:
#
#   1. scan_dat_candidates (checkpoint) — scans the whole .dat file directly
#      and nominates every period bin meeting the candidate_scan filters.
#      No manual pause/curation step — every passing bin proceeds.
#   2. extract_by_period — per bin, pull consensus motifs (TRF .dat field 14)
#      for entries within period +/- window.
#   3. rank_family_clusters — per bin, iteratively peel off coherent
#      sequence clusters (rotation/strand-aware anchor matching, same trick
#      as before) and rank them by size: O(rounds x N) per bin, not O(N^2)
#      all-vs-all. Produces one consensus per cluster (rank 1 = most
#      support, rank 2, ...). Sequences that never join a cluster of
#      min_cluster_size are dropped, not reported.
#   4. build_all_clusters_table — aggregates every bin's ranked clusters
#      into one master table.
#   5. cross_cluster_comparison — pairwise rotation/strand-aware comparison
#      across EVERY cluster found (all bins, all ranks — including within
#      the same bin), flagging near-identical or integer-multiple pairs.
#   6. resolve_redundancy — processes clusters in descending support order;
#      each cluster is checked only against already-confirmed winners (never
#      against another loser, never chained transitively through an
#      intermediate). If flagged against a winner it's marked redundant with
#      that winner specifically; otherwise it becomes a winner itself.
#      Nothing is dropped — losers are kept in the summary table and the
#      RepeatMasker library, marked redundant (both in the summary table and
#      in the library's classification field), rather than removed.
#   7. known-repeat screening (prepare_known_repeat_query ->
#      run_repeatmasker_known_screen -> parse_known_repeat_hits) — runs
#      RepeatMasker -species against OUR OWN discovered consensus motifs
#      (not the assembly) to check whether any match a previously
#      characterized repeat family in Dfam/RepBase. A targeted, much
#      cheaper "is this already known" check than screening the whole
#      assembly. Uses a fresh, pipeline-owned RepeatMasker install (conda
#      only, see workflow/envs/repeatmasker.yaml — floor-pinned to the
#      latest release. The FamDB library it needs does NOT come bundled or
#      auto-downloaded on this version chain, so the first run after a
#      fresh env build will fail with "FamDB data directory not found" —
#      this is expected, not a bug; see the README's "Known-repeat
#      screening" section for the one-time manual setup), independent of
#      whatever version/library may already be on the cluster.
#
# This requires Snakemake checkpoints (checkpoint scan_dat_candidates):
# the bin list isn't known until the scan actually runs, so downstream rules
# use input functions that call checkpoints.scan_dat_candidates.get(...) to
# discover it at DAG-re-evaluation time, rather than a static manifest.
#
# Usage: snakemake -s Snakefile --configfile config/config.yaml --cores <N> \
#          --use-envmodules --retries 3
# See README.md for full documentation and cluster submission instructions.
#
# Author:  KM
# Created: 2026-08
# =============================================================================

wildcard_constraints:
    motif="[^/]+",


def get_manifest(wildcards):
    """Parse the scan checkpoint's manifest into {name: (period, window)}.
    Every rule/function that needs the bin list or a bin's period/window
    calls this, which forces the checkpoint to complete first."""
    path = checkpoints.scan_dat_candidates.get(**wildcards).output.manifest
    entries = {}
    with open(path) as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                fields = line.rstrip("\n").split("\t")
                entries[fields[0]] = (int(fields[1]), int(fields[2]))
    return entries


def all_motifs(wildcards):
    return list(get_manifest(wildcards).keys())


def motif_period(wildcards):
    return get_manifest(wildcards)[wildcards.motif][0]


def motif_window(wildcards):
    return get_manifest(wildcards)[wildcards.motif][1]


def all_motif_periods(wildcards):
    m = get_manifest(wildcards)
    return [m[name][0] for name in m]


def all_motif_windows(wildcards):
    m = get_manifest(wildcards)
    return [m[name][1] for name in m]


rule all:
    input:
        "results/repeatmasker_custom_lib.fasta",
        "results/summary_table.tsv",
        "results/known_repeat_hits.tsv"


checkpoint scan_dat_candidates:
    # Scans the whole .dat file directly; the bin list downstream rules use
    # is determined here, at runtime, not known when the DAG is built.
    input:
        dat=config["trf_dat"]
    output:
        tsv="results/candidate_periods.tsv",
        raw_tsv="results/candidate_periods_raw.tsv",
        manifest="results/candidate_periods_manifest.tsv"
    log:
        "results/logs/scan_dat_candidates.log"
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
        dat=config["trf_dat"],
        manifest=lambda wc: checkpoints.scan_dat_candidates.get(**wc).output.manifest
    output:
        fasta="results/{motif}/01_raw_consensus.fasta",
        tsv="results/{motif}/01_raw_consensus.tsv"
    params:
        period=motif_period,
        window=motif_window
    log:
        "results/logs/{motif}/extract_by_period.log"
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
        fasta="results/{motif}/01_raw_consensus.fasta"
    output:
        clusters_dir=directory("results/{motif}/ranked_clusters"),
        summary="results/{motif}/ranked_clusters_summary.tsv"
    params:
        motif_name=lambda wc: wc.motif
    log:
        "results/logs/{motif}/rank_family_clusters.log"
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
        "--in-fasta {input.fasta} --out-dir {output.clusters_dir} "
        "--out-summary {output.summary} --motif-name {params.motif_name} "
        "--min-identity {config[cluster_ranking][min_cluster_identity]} "
        "--max-gap-fraction {config[cluster_ranking][max_gap_fraction]} "
        "--min-cluster-size {config[cluster_ranking][min_cluster_size]} "
        "--max-rounds {config[cluster_ranking][max_rounds]} > {log} 2>&1"


rule build_all_clusters_table:
    input:
        summaries=lambda wc: expand("results/{motif}/ranked_clusters_summary.tsv", motif=all_motifs(wc)),
        dirs=lambda wc: expand("results/{motif}/ranked_clusters", motif=all_motifs(wc))
    output:
        "results/all_clusters.tsv"
    params:
        motifs=lambda wc: all_motifs(wc),
        periods=lambda wc: all_motif_periods(wc),
        windows=lambda wc: all_motif_windows(wc)
    log:
        "results/logs/build_all_clusters_table.log"
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
    # across every bin — safe at runtime since Snakemake guarantees the
    # `directory()` outputs from rank_family_clusters are fully materialized
    # before this rule executes (both are declared inputs).
    input:
        summaries=lambda wc: expand("results/{motif}/ranked_clusters_summary.tsv", motif=all_motifs(wc)),
        dirs=lambda wc: expand("results/{motif}/ranked_clusters", motif=all_motifs(wc))
    output:
        "results/all_clusters_consensus.fasta"
    log:
        "results/logs/combine_cluster_fastas.log"
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


rule cross_cluster_comparison:
    input:
        fasta="results/all_clusters_consensus.fasta"
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
        clusters_tsv="results/all_clusters.tsv",
        comparison_tsv="results/cross_cluster_comparison.tsv",
        consensus_fasta="results/all_clusters_consensus.fasta"
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


# --- Step 7: known-repeat screening ---
# Screens our own discovered consensus motifs (not the assembly) against a
# Dfam/RepBase library via RepeatMasker -species, to check whether any
# match a previously characterized repeat family.

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
    # Deliberately conda-only, no envmodules: fallback — unlike every other
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
        outdir="results/known_repeat_screen"
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
        "RepeatMasker -species '{params.species}' -pa {threads} -uncurated "
        "-dir {params.outdir} {input} > {log} 2>&1 && touch {output}"


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
