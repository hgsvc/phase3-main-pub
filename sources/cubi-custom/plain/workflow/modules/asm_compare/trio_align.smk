

EVAL_ALIGN_PRECISION = ["exact", "struct-50"]


rule align_verkko_parent_to_child:
    input:
        child = WORKDIR_EVAL.joinpath(
            "proc/05-preprocess/merge_tag_asm",
            "{child}.vrk-ps-sseq.asm-mrg-tag.fasta"
        ),
        chl_idx = WORKDIR_EVAL.joinpath(
            "proc/05-preprocess/merge_tag_asm",
            "{child}.vrk-ps-sseq.asm-mrg-tag.fasta.fai"
        ),
        parent = WORKDIR_EVAL.joinpath(
            "proc/05-preprocess/merge_tag_asm",
            "{parent}.vrk-ps-sseq.asm-mrg-tag.fasta"
        ),
        par_idx = WORKDIR_EVAL.joinpath(
            "proc/05-preprocess/merge_tag_asm",
            "{parent}.vrk-ps-sseq.asm-mrg-tag.fasta.fai"
        ),
    output:
        paf = DIR_PROC.joinpath(
            "asm_compare", "trio_align", "raw",
            "{parent}-to-{child}.vrk-ps-sseq.paf.gz"
        )
    conda:
        DIR_ENVS.joinpath("minimap.yaml")
    threads: 12
    resources:
        mem_mb=lambda wildcards, attempt: 65536 * attempt,
        time_hrs=lambda wildcards, attempt: 4 * attempt
    shell:
        "minimap2 -t {threads} -x asm5 --eqx --MD -c --secondary=no "
        "{input.child} {input.parent} | pigz > {output.paf}"


rule normalize_trio_align_parent_to_child:
    input:
        paf = rules.align_verkko_parent_to_child.output.paf
    output:
        tsv = DIR_PROC.joinpath(
            "asm_compare", "trio_align", "norm",
            "{parent}-to-{child}.vrk-ps-sseq.norm-paf.tsv.gz"
        )
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    params:
        script=DIR_SCRIPTS.joinpath("normalize_paf.py").resolve(strict=True)
    shell:
        "{params.script} --input {input.paf} --output {output.tsv}"


rule evaluate_normalized_trio_pafs:
    input:
        paf = DIR_PROC.joinpath(
            "asm_compare", "trio_align", "norm",
            "{parent}-to-{child}.vrk-ps-sseq.norm-paf.tsv.gz"
        ),
        gsize = select_genome_size_file
    output:
        stats = DIR_RES.joinpath(
            "asm_compare", "trio",
            "statistics",
            "{parent}-to-{child}.vrk-ps-sseq.mapq-{min_mapq}.seq-{min_seq_len}.aln-{min_aln_len}.{precision}.stats.tsv"
        ),
        aln_regions = DIR_RES.joinpath(
            "asm_compare", "trio",
            "regions",
            "{parent}-to-{child}.vrk-ps-sseq.mapq-{min_mapq}.seq-{min_seq_len}.aln-{min_aln_len}.{precision}.aligned-regions.tsv.gz"
        ),
        sup_regions = DIR_RES.joinpath(
            "asm_compare", "trio",
            "regions",
            "{parent}-to-{child}.vrk-ps-sseq.mapq-{min_mapq}.seq-{min_seq_len}.aln-{min_aln_len}.{precision}.support-regions.tsv.gz"
        ),
        mrg_regions = DIR_RES.joinpath(
            "asm_compare", "trio",
            "regions",
            "{parent}-to-{child}.vrk-ps-sseq.mapq-{min_mapq}.seq-{min_seq_len}.aln-{min_aln_len}.{precision}.merged-regions.tsv.gz"
        )
    log:
        DIR_RES.joinpath(
            "asm_compare", "trio",
            "statistics",
            "{parent}-to-{child}.vrk-ps-sseq.mapq-{min_mapq}.seq-{min_seq_len}.aln-{min_aln_len}.{precision}.aln-filter.log"
        )
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt * attempt
    params:
        script=DIR_SCRIPTS.joinpath("asm_compare", "eval_paf.py").resolve(strict=True),
        min_aln_len=lambda wildcards: to_int(wildcards.min_aln_len),
        min_seq_len=lambda wildcards: to_int(wildcards.min_seq_len),
        precision=lambda wildcards: parse_precision_param(wildcards.precision, "setting"),
        struct_t=lambda wildcards: parse_precision_param(wildcards.precision, "threshold"),
    shell:
        "{params.script} --verbose --paf-file {input.paf} "
        "--genome-size {input.gsize} "
        "--precision {params.precision} "
        "--retain-mapq-geq {wildcards.min_mapq} "
        "--retain-seqlen-geq {params.min_seq_len} "
        "--retain-alnlen-geq {params.min_aln_len} "
        "--struct-err-geq {params.struct_t} "
        "--split-seq-tag '.' "
        "--out-seq-stats {output.stats} "
        "--out-aligned-regions {output.aln_regions} "
        "--out-support-regions {output.sup_regions} "
        "--out-merged-regions {output.mrg_regions} "
        "&> {log}"


rule postprocess_merged_trio_regions:
    input:
        mrg_regions = rules.evaluate_normalized_trio_pafs.output.mrg_regions,
        gsize = select_genome_size_file
    output:
        mrg_stats = DIR_RES.joinpath(
            "asm_compare", "trio",
            "statistics",
            "{parent}-to-{child}.vrk-ps-sseq.mapq-{min_mapq}.seq-{min_seq_len}.aln-{min_aln_len}.{precision}.mrg-regions-stats.tsv"
        )
    run:
        import pandas as pd

        seq_lens = dict()
        with open(input.gsize, "r") as listing:
            for line in listing:
                seq_name, length = line.split()[:2]
                seq_lens[seq_name] = int(length)
                # genome size file is fasta index by default;
                # that contains tagged sequence names
                seq_name, _ = seq_name.rsplit(".", 1)
                seq_lens[seq_name] = int(length)

        df = pd.read_csv(input.mrg_regions, sep="\t", header=0)
        df["length"] = df["end"] - df["start"]

        pct_stats = df["length"].describe(percentiles=[0.25,0.5,0.75,0.95])
        pct_stats.rename(
            {
                "count": "num_regions", "mean": "length_mean",
                "std": "stddev", "min": "length_pct_0",
                "25%": "length_pct_25", "50%": "length_pct_50", "75%": "length_pct_75",
                "95%": "length_pct_95", "max": "length_pct_100"
            },
            inplace=True
        )
        num_seqs = df["seq_name"].nunique()
        known_seqs = set(df["seq_name"].unique())
        total_length = sum(seq_lens[sn] for sn in known_seqs)
        support_length = df["length"].sum()
        support_pct = round(support_length/total_length * 100, 2)

        stats_row = [
            wildcards.parent, wildcards.child, wildcards.precision,
            num_seqs, total_length, support_length, support_pct,
            int(pct_stats["num_regions"]), int(round(pct_stats["length_mean"], 0)),
            int(round(pct_stats["length_pct_25"],0)), int(round(pct_stats["length_pct_50"],0)),
            int(round(pct_stats["length_pct_75"],0)), int(round(pct_stats["length_pct_95"],0)),
            int(round(pct_stats["length_pct_100"],0))
        ]

        header = [
            "parent", "child", "precision",
            "num_contigs", "total_length", "supported_length", "supported_pct",
            "num_support_regions", "region_length_mean", "region_length_pct25",
            "region_length_pct50", "region_length_pct75", "region_length_pct95",
            "region_length_pct100"
        ]

        with open(output.mrg_stats, "w") as dump:
            _ = dump.write("\t".join(header) + "\n")
            _ = dump.write("\t".join(map(str, stats_row)) + "\n")
    # END OF RUN BLOCK


TRIO_MAP = {
    "HG00514": ["HG00512", "HG00513"],
    "HG00733": ["HG00731", "HG00732"],
    "NA19240": ["NA19238", "NA19239"]
}

FAMILY_MAP = {
    "HG00514": "SH032-CHS",
    "HG00733": "PR05-PUR",
    "NA19240": "Y117-YRI"
}


rule concat_merged_trio_region_stats:
    input:
        tables = lambda wildcards: expand(
            rules.postprocess_merged_trio_regions.output.mrg_stats,
            parent=TRIO_MAP[wildcards.child],
            precision=EVAL_ALIGN_PRECISION,
            allow_missing=True
        )
    output:
        summary = DIR_RES.joinpath(
            "asm_compare", "trio",
            "statistics", "summary",
            "{child}.vrk-ps-sseq.mapq-{min_mapq}.seq-{min_seq_len}.aln-{min_aln_len}.merged-regions.tsv"
        )
    run:
        import pandas as pd
        import pathlib as pl

        concat = []
        for table_file in input.tables:
            df = pd.read_csv(table_file, header=0, sep="\t")
            concat.append(df)
        concat = pd.concat(concat, axis=0, ignore_index=False)
        concat.sort_values(["parent", "precision"], inplace=True)

        concat.to_csv(output.summary, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


localrules: create_parental_summary
rule create_parental_summary:
    input:
        parents = lambda wildcards: expand(
            rules.evaluate_normalized_trio_pafs.output.stats,
            parent=TRIO_MAP[wildcards.child],
            allow_missing=True
        )
    output:
        table = DIR_RES.joinpath(
            "asm_compare", "trio",
            "statistics", "summary",
            "{child}.vrk-ps-sseq.mapq-{min_mapq}.seq-{min_seq_len}.aln-{min_aln_len}.{precision}.hap-support-{size_cutoff}.tsv"
        ),
        min_summ = DIR_RES.joinpath(
            "asm_compare", "trio",
            "statistics", "summary",
            "{child}.vrk-ps-sseq.mapq-{min_mapq}.seq-{min_seq_len}.aln-{min_aln_len}.{precision}.hap-support-{size_cutoff}.summary.tsv"
        )
    params:
        size_cutoff = lambda wildcards: to_int(wildcards.size_cutoff)
    run:
        import pandas as pd
        import pathlib as pl

        size_lower_bound = params.size_cutoff
        merged = None
        parent1 = None
        parent2 = None
        for table_file in sorted(input.parents):
            parent = pl.Path(table_file).name.split("-")[0]
            df = pd.read_csv(table_file, sep="\t", header=0)
            df = df.loc[df["seq_size"] > size_lower_bound, :].copy()
            df = df[["seq_name", "seq_tag", "seq_size", "merged_sup_length", "merged_sup_pct"]].copy()
            df.rename(
                {
                    "merged_sup_length": f"{parent}_support_bp", "merged_sup_pct": f"{parent}_support_pct"
                }, axis=1, inplace=True
            )
            if parent1 is None:
                parent1 = parent
                p1_column = f"{parent}_support_pct"
            else:
                parent2 = parent
                p2_column = f"{parent}_support_pct"
            df["child"] = wildcards.child
            if merged is None:
                merged = df.copy()
            else:
                merged = merged.merge(df, on=["child", "seq_name", "seq_tag", "seq_size"], how="outer")
        merged.sort_values(["seq_name", "seq_size"], inplace=True)

        merged["parent_haplotype"] = "undetermined"
        is_unassigned = merged["seq_tag"] == "unassigned"
        if is_unassigned.any():
            merged.loc[is_unassigned, "parent_haplotype"] = "unphased"
        merged.loc[merged[p1_column] > merged[p2_column], "parent_haplotype"] = parent1
        merged.loc[merged[p1_column] < merged[p2_column], "parent_haplotype"] = parent2

        merged = merged[
            [
                "child", "seq_name", "seq_tag", "seq_size", "parent_haplotype",
                f"{parent1}_support_bp", f"{parent1}_support_pct",
                f"{parent2}_support_bp", f"{parent2}_support_pct"
            ]
        ].copy()
        merged.to_csv(output.table, sep="\t", header=True, index=False)

        minimal_summary = merged.groupby("parent_haplotype").agg(
            haplotype_seq=pd.NamedAgg(column="seq_name", aggfunc="nunique"),
            haplotype_length=pd.NamedAgg(column="seq_size", aggfunc="sum"),
            parent1_pct_support_median=pd.NamedAgg(column=p1_column, aggfunc="median"),
            parent1_pct_support_min=pd.NamedAgg(column=p1_column, aggfunc="min"),
            parent1_pct_support_max=pd.NamedAgg(column=p1_column, aggfunc="max"),
            parent2_pct_support_median=pd.NamedAgg(column=p2_column, aggfunc="median"),
            parent2_pct_support_min=pd.NamedAgg(column=p2_column, aggfunc="min"),
            parent2_pct_support_max=pd.NamedAgg(column=p2_column, aggfunc="max"),
        )
        minimal_summary.insert(0, "child", wildcards.child)
        minimal_summary.to_csv(output.min_summ, sep="\t", header=True, index=True, index_label="parent_haplotype")
    # END OF RUN BLOCK


localrules: merge_parental_summaries
rule merge_parental_summaries:
    input:
        tables = lambda wildcards: expand(
            rules.create_parental_summary.output.table,
            child=list(TRIO_MAP.keys()),
            allow_missing=True
        ),
        minimals = lambda wildcards: expand(
            rules.create_parental_summary.output.min_summ,
            child=list(TRIO_MAP.keys()),
            allow_missing=True
        )
    output:
        minimal = DIR_RES.joinpath(
            "asm_compare", "trio",
            "statistics", "merged",
            "TRIOS.vrk-ps-sseq.mapq-{min_mapq}.seq-{min_seq_len}.aln-{min_aln_len}.{precision}.hap-support-{size_cutoff}.summary.tsv"
        ),
        reported_number = DIR_RES.joinpath(
            "asm_compare", "trio",
            "statistics", "merged",
            "TRIOS.vrk-ps-sseq.mapq-{min_mapq}.seq-{min_seq_len}.aln-{min_aln_len}.{precision}.hap-support-{size_cutoff}.median.txt"
        )
    run:
        import pandas as pd
        import statistics
        import pathlib as pl
        import io

        support_values = []
        family_medians = io.StringIO()
        for table in input.tables:
            family_values = []
            child = pl.Path(table).name.split(".")[0]
            family = FAMILY_MAP[child]
            df = pd.read_csv(table, sep="\t", header=0)
            for parent, stats in df.groupby("parent_haplotype"):
                try:
                    values = stats[f"{parent}_support_pct"].values
                    support_values.extend(list(values))
                    family_values.extend(list(values))
                except KeyError:
                    support_values.extend([0] * stats.shape[0])
                    family_values.extend([0] * stats.shape[0])
            fam_med = round(statistics.median(family_values), 2)
            family_medians.write(f"{family}_median_parental_support\t{fam_med}\n")
        median = round(statistics.median(support_values), 2)
        with open(output.reported_number, "w") as dump:
            dump.write(family_medians.getvalue())
            dump.write(f"trio_median_parental_support\t{median}\n")

        concat = []
        for summary in input.minimals:
            df = pd.read_csv(summary, sep="\t", header=0)
            concat.append(df)

        concat = pd.concat(concat, axis=0, ignore_index=False)
        concat.sort_values(["child", "parent_haplotype"], inplace=True)

        concat.to_csv(output.minimal, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule run_all_asm_trio_align:
    input:
        summary = expand(
            rules.concat_merged_trio_region_stats.output.summary,
            child=["HG00733", "HG00514", "NA19240"],
            min_mapq=[1],
            min_seq_len=["0", "100k"],
            min_aln_len=["0", "10k"]
        ),
        min_summ = expand(
            rules.create_parental_summary.output.table,
            child=["HG00733", "HG00514", "NA19240"],
            min_mapq=[1],
            min_seq_len=["100k"],
            min_aln_len=["10k"],
            precision=EVAL_ALIGN_PRECISION,
            size_cutoff=["1M"]
        ),
        number = expand(
            rules.merge_parental_summaries.output.reported_number,
            min_mapq=[1],
            min_seq_len=["100k"],
            min_aln_len=["10k"],
            precision=EVAL_ALIGN_PRECISION,
            size_cutoff=["1M"]
        )
