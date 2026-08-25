#!/usr/bin/env python
"""
filter_and_pool_clusters.py

Stage 3 of the multi-entry pipeline: reads every manifest entry's
all_clusters.tsv (build_all_clusters_table.py's per-entry output, already
aggregated across that entry's period bins), drops clusters whose
total_copy_number falls below --min-total-copy-number, and pools every
surviving cluster from every entry into one table and one concatenated
consensus FASTA.

Every surviving cluster is renamed cluster_uid = "{entry_id}__{label}" so
labels that collide across entries (e.g. two entries both producing a
"cand171_rank1_n..." label) stay distinct once pooled. gc_content is
carried through from all_clusters.tsv here specifically so downstream
per-family provenance (resolve_redundancy.py) has it without needing to
re-open every entry's raw table.

--clusters-tsv/--fastas/--entry-ids/--individuals/--assembly-methods/
--phasing-statuses are all positionally paired, one entry per index (same
convention build_all_clusters_table.py uses for --motifs/--periods/--windows).
"""
import argparse
import sys


def read_tsv(path):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = [dict(zip(header, line.rstrip("\n").split("\t"))) for line in f if line.strip()]
    return rows


def read_multi_fasta(path):
    """Label is each record's header up to the first space — matches
    rank_family_clusters.py's own FASTA header convention."""
    seqs = {}
    header, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(seq)
                header = line[1:].split()[0] if line[1:].split() else line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if header is not None:
            seqs[header] = "".join(seq)
    return seqs


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--clusters-tsv", nargs="+", required=True,
                     help="Each entry's all_clusters.tsv (pre-filter), one per entry")
    ap.add_argument("--fastas", nargs="+", required=True,
                     help="Each entry's all_clusters_consensus.fasta, positionally paired "
                          "with --clusters-tsv")
    ap.add_argument("--entry-ids", nargs="+", required=True)
    ap.add_argument("--individuals", nargs="+", required=True)
    ap.add_argument("--assembly-methods", nargs="+", required=True)
    ap.add_argument("--phasing-statuses", nargs="+", required=True)
    ap.add_argument("--min-total-copy-number", type=float, required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-fasta", required=True)
    args = ap.parse_args()

    n = len(args.entry_ids)
    for name, lst in [
        ("clusters-tsv", args.clusters_tsv), ("fastas", args.fastas),
        ("individuals", args.individuals), ("assembly-methods", args.assembly_methods),
        ("phasing-statuses", args.phasing_statuses),
    ]:
        if len(lst) != n:
            sys.exit(f"ERROR: --{name} has {len(lst)} entries, expected {n} (one per --entry-ids)")

    out_columns = ["cluster_uid", "entry_id", "individual", "assembly_method", "phasing_status",
                   "source_cluster_label", "consensus_length", "gc_content",
                   "n_input_sequences", "total_copy_number"]

    rows_out = []
    fasta_records = []
    n_seen = 0
    n_kept = 0

    for i in range(n):
        entry_id = args.entry_ids[i]
        individual = args.individuals[i]
        assembly_method = args.assembly_methods[i]
        phasing_status = args.phasing_statuses[i]
        seqs = read_multi_fasta(args.fastas[i])

        for row in read_tsv(args.clusters_tsv[i]):
            n_seen += 1
            total_copy_number = float(row["total_copy_number"])
            if total_copy_number < args.min_total_copy_number:
                continue
            n_kept += 1

            label = row["label"]
            cluster_uid = f"{entry_id}__{label}"
            seq = seqs.get(label)
            if seq is None:
                sys.exit(
                    f"ERROR: cluster '{label}' from {args.clusters_tsv[i]} has no matching "
                    f"record in {args.fastas[i]} (entry_id={entry_id}) — clusters-tsv/fastas "
                    f"pair is inconsistent for this entry"
                )

            rows_out.append({
                "cluster_uid": cluster_uid, "entry_id": entry_id, "individual": individual,
                "assembly_method": assembly_method, "phasing_status": phasing_status,
                "source_cluster_label": label, "consensus_length": row["consensus_length"],
                "gc_content": row["gc_content"], "n_input_sequences": row["n_input_sequences"],
                "total_copy_number": row["total_copy_number"],
            })
            fasta_records.append((cluster_uid, row["n_input_sequences"], row["total_copy_number"], seq))

    with open(args.out_tsv, "w") as out:
        out.write("\t".join(out_columns) + "\n")
        for row in rows_out:
            out.write("\t".join(str(row[c]) for c in out_columns) + "\n")

    with open(args.out_fasta, "w") as out:
        for cluster_uid, n_input_sequences, total_copy_number, seq in fasta_records:
            out.write(
                f">{cluster_uid} n_input_sequences={n_input_sequences} "
                f"total_copy_number={total_copy_number}\n{seq}\n"
            )

    sys.stderr.write(
        f"[filter_and_pool_clusters] {n} entries, {n_seen} clusters seen, {n_kept} survived "
        f"total_copy_number >= {args.min_total_copy_number}, {n_seen - n_kept} filtered out\n"
    )


if __name__ == "__main__":
    main()
