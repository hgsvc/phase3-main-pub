
rule segdups_to_bed:
    input:
        sd95 = rules.split_segdup_annotation.output.sd95,
        sd98 = rules.split_segdup_annotation.output.sd98
    output:
        bed_95 = temp(DIR_PROC.joinpath(
            "regions", "gaps", "label_segdup",
            "{sample}.sd-095.bed"
        )),
        bed_98 = temp(DIR_PROC.joinpath(
            "regions", "gaps", "label_segdup",
            "{sample}.sd-098.bed"
        ))
    run:
        import pandas as pd
        import hashlib as hl

        def compute_row_id(row):

            row_hash = hl.md5(f"{row.seq}{row.start}{row.end}".encode("utf-8")).hexdigest()
            row_name = f"{row_hash}|{row.label}"
            return row_name


        for input_file, output_file, region_label in zip(input, output, ["SD95", "SD98"]):
            df = pd.read_csv(input_file, sep="\t", header=0)
            df["label"] = region_label
            df["name"] = df.apply(compute_row_id, axis=1)
            df.rename({"seq": "#contig"}, axis=1, inplace=True)
            df.sort_values(["#contig", "start", "end"], inplace=True)
            df.to_csv(output_file, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule merge_overlapping_segdups:
    input:
        bed = DIR_PROC.joinpath(
            "regions", "gaps", "label_segdup",
            "{sample}.sd-{pct_id}.bed"
        ),
    output:
        bed = DIR_PROC.joinpath(
            "regions", "gaps", "label_segdup",
            "{sample}.sd-{pct_id}.mrg.bed"
        ),
    wildcard_constraints:
        pct_id = "(095|098)"
    conda:
        DIR_ENVS.joinpath("bedtools.yaml")
    shell:
        "bedtools merge -c 5 -o first -i {input.bed} > {output.bed}"


rule annotate_gaps_with_segdups:
    """Annotate assemblies (assembly coordinate space)
    with segdup annotation. The input file in ref coordinates
    is just here for enforcing.
    """
    input:
        qry_view = WORKDIR_EVAL.joinpath(
            "results", "regions", "{sample}",
            "{sample}.asm-{asm_unit}.{refgenome}.ctg-aln-gap.asm-coord.bed"
        ),
        trg_view = WORKDIR_EVAL.joinpath(
            "results", "regions", "{sample}",
            "{sample}.asm-{asm_unit}.{refgenome}.ctg-aln-gap.ref-coord.bed"
        ),
        segdups = rules.merge_overlapping_segdups.output.bed
    output:
        isect = DIR_PROC.joinpath(
            "regions", "gaps", "intersections",
            "{sample}",
            "{sample}.asm-{asm_unit}.{refgenome}.isect-sd{pct_id}.tsv"
        )
    conda:
        DIR_ENVS.joinpath("bedtools.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    shell:
        "bedtools intersect -wao -a {input.qry_view} -b {input.segdups} > {output.isect}"


rule simplify_segdup_gap_intersection:
    input:
        table = rules.annotate_gaps_with_segdups.output.isect
    output:
        tsv = DIR_RES.joinpath(
            "regions", "gaps", "norm_tables",
            "{sample}",
            "{sample}.asm-{asm_unit}.{refgenome}.asm-coord.segdup-{pct_id}.tsv"
        )
    run:
        import pandas as pd
        header = ["aln_seq", "aln_start", "aln_end", "aln_label", "aln_length", "aln_info"]
        header += ["aln_base_block", "aln_coarse_block"]
        header += ["sd_seq", "sd_start", "sd_end", "sd_id"]
        header += ["overlap_bp"]

        df = pd.read_csv(input.table, sep="\t", header=None, names=header)
        df = df.loc[df["overlap_bp"] > 0, :].copy()
        df["sd_length"] = df["sd_end"] - df["sd_start"]
        df["overlap_pct"] = (df["overlap_bp"] / df["sd_length"] * 100).round(2)
        # drop all SDs that are inside an aligned block - not interesting
        is_aligned = df["aln_label"] == "ALN"
        is_complete = df["overlap_pct"] > 99.9
        select = is_aligned & is_complete
        df = df.loc[~select, :].copy()
        df.drop("aln_length", axis=1, inplace=True)
        df["sample"] = wildcards.sample
        df["asm_unit"] = wildcards.asm_unit
        # the following does not work unfortunately
        #assert "NO-BLOCK-INFO" not in set(df["aln_base_block"].unique())
        df.sort_values(["aln_seq", "aln_start", "aln_end"], inplace=True)
        df.to_csv(output.tsv, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


localrules: hprc_gaps_to_bed
rule hprc_gaps_to_bed:
    input:
        tsv = HPRC_COMMON_GAPS
    output:
        bed = DIR_PROC.joinpath(
            "regions", "gaps", "norm_tables",
            "hprc_common_gaps.bed"
        )
    run:
        import pandas as pd
        df = pd.read_csv(input.tsv, sep="\t", header=0)
        df.sort_values(["chrom", "start", "end"], inplace=True)
        if (df["start"] == df["end"]).any():
            # assume 1-based coordinates
            df["start"] -= 1
            assert (df["start"] >= 0).all()
        df.rename({"chrom": "#chrom"}, axis=1, inplace=True)
        df.to_csv(output.bed, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule annotate_gaps_with_hprc_gaps:
    """NB: HPRC gaps only exist in T2T space
    """
    input:
        qry_view = WORKDIR_EVAL.joinpath(
            "results", "regions", "{sample}",
            "{sample}.asm-{asm_unit}.{refgenome}.ctg-aln-gap.asm-coord.bed"
        ),
        trg_view = WORKDIR_EVAL.joinpath(
            "results", "regions", "{sample}",
            "{sample}.asm-{asm_unit}.{refgenome}.ctg-aln-gap.ref-coord.bed"
        ),
        gaps = rules.hprc_gaps_to_bed.output.bed
    output:
        isect = DIR_PROC.joinpath(
            "regions", "gaps", "intersections",
            "{sample}",
            "{sample}.asm-{asm_unit}.{refgenome}.isect-hprcgaps.tsv"
        )
    conda:
        DIR_ENVS.joinpath("bedtools.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    shell:
        "bedtools intersect -wao -a {input.trg_view} -b {input.gaps} > {output.isect}"


rule simplify_hprc_gap_intersection:
    input:
        table = rules.annotate_gaps_with_hprc_gaps.output.isect
    output:
        tsv = DIR_RES.joinpath(
            "regions", "gaps", "norm_tables",
            "{sample}",
            "{sample}.asm-{asm_unit}.{refgenome}.ref-coord.hprc-gaps.tsv"
        )
    run:
        import pandas as pd
        header = ["aln_seq", "aln_start", "aln_end", "aln_label", "aln_length", "aln_info"]
        header += ["aln_base_block", "aln_coarse_block"]
        header += ["gap_chrom", "gap_start", "gap_end", "gap_id", "hprc_haps"]
        header += ["gap_sd_assoc", "gap_cov_drop", "gap_win_start", "gap_win_end", "overlap_bp"]

        df = pd.read_csv(input.table, sep="\t", header=None, names=header)
        df = df.loc[df["overlap_bp"] > 0, :].copy()
        df.drop(["hprc_haps", "gap_sd_assoc", "gap_cov_drop", "gap_win_start", "gap_win_end", "aln_length"], axis=1, inplace=True)
        # the following does not work unfortunately
        #assert "NO-BLOCK-INFO" not in set(df["aln_base_block"].unique())
        df["sample"] = wildcards.sample
        df["asm_unit"] = wildcards.asm_unit
        df["gap_length"] = df["gap_end"] - df["gap_start"]
        df["overlap_pct"] = (df["overlap_bp"] / df["gap_length"] * 100).round(2)
        df.sort_values(["aln_seq", "aln_start", "aln_end"], inplace=True)
        df.to_csv(output.tsv, sep="\t", header=True, index=False)
    # END OR RUN BLOCK


rule filter_nucfreq_regions_to_bed:
    input:
        tsv = WORKDIR_EVAL.joinpath(
            "results", "regions",
            "{sample}",
            "{sample}.nucfreq.covann.tsv.gz"
        )
    output:
        hap1 = DIR_PROC.joinpath(
            "regions", "gaps", "norm_tables",
            "{sample}",
            "{sample}.asm-hap1.nucfreq.covann-filtered.bed"
        ),
        hap2 = DIR_PROC.joinpath(
            "regions", "gaps", "norm_tables",
            "{sample}",
            "{sample}.asm-hap2.nucfreq.covann-filtered.bed"
        ),
    run:
        import pandas as pd

        df = pd.read_csv(input.tsv, sep="\t", header=0)
        threshold = 125
        df = df.loc[df["hifi_pct_median_cov"] > threshold, :].copy()
        df.rename({"contig": "#contig"}, axis=1, inplace=True)

        hap1 = df.loc[df["asm_unit"] == "hap1", :].copy()
        hap2 = df.loc[df["asm_unit"] == "hap2", :].copy()

        hap1.to_csv(output.hap1, sep="\t", header=True, index=False)
        hap2.to_csv(output.hap2, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule split_normalized_flagger:
    input:
        tsv = rules.normalize_flagger_annotation.output.bed_like
    output:
        hap1 = DIR_PROC.joinpath(
            "regions", "gaps", "norm_tables",
            "{sample}",
            "{sample}.asm-hap1.flagger-regions.bed"
        ),
        hap2 = DIR_PROC.joinpath(
            "regions", "gaps", "norm_tables",
            "{sample}",
            "{sample}.asm-hap2.flagger-regions.bed"
        ),
        un = DIR_PROC.joinpath(
            "regions", "gaps", "norm_tables",
            "{sample}",
            "{sample}.asm-unassigned.flagger-regions.bed"
        ),
    run:
        import pandas as pd

        df = pd.read_csv(input.tsv, sep="\t", header=0)

        def assign_unit(contig_name):
            if "haplotype1" in contig_name:
                return "hap1"
            elif "haplotype2" in contig_name:
                return "hap2"
            elif "unassigned" in contig_name:
                return "unassigned"
            else:
                raise

        df["asm_unit"] = df["chrom"].apply(assign_unit)
        # for gap / QC annotation, keep only labels that
        # indicate any type of error
        df = df.loc[df["name"] != "Hap", :].copy()
        df.rename({"chrom": "#contig"}, axis=1, inplace=True)
        hap1 = df.loc[df["asm_unit"] == "hap1", :].copy()
        hap1.to_csv(output.hap1, sep="\t", header=True, index=False)

        hap2 = df.loc[df["asm_unit"] == "hap2", :].copy()
        hap2.to_csv(output.hap2, sep="\t", header=True, index=False)

        un = df.loc[df["asm_unit"] == "unassigned", :].copy()
        un.to_csv(output.un, sep="\t", header=True, index=False)
    # END OR RUN BLOCK


rule annotate_gaps_with_qc_info:
    input:
        gaps = rules.hprc_gaps_to_bed.output.bed,
        isect_table = rules.simplify_hprc_gap_intersection.output.tsv,
        flagger = DIR_PROC.joinpath(
            "regions", "gaps", "norm_tables",
            "{sample}",
            "{sample}.asm-{asm_unit}.flagger-regions.bed"
        ),
        nucfreq = DIR_PROC.joinpath(
            "regions", "gaps", "norm_tables",
            "{sample}",
            "{sample}.asm-{asm_unit}.nucfreq.covann-filtered.bed"
        ),
        karyotypes = WORKDIR_EVAL.joinpath(
            "results/reports/ref_chrom_assign",
            "karyo-est.hgsvc3-verkko.tsv"
        ),
        chrom_assign = WORKDIR_EVAL.joinpath(
            "results/reports/ref_chrom_assign",
            "{sample}.asm-{asm_unit}.{refgenome}.chrom-assign-by-query.tsv"
        )
    output:
        table = DIR_RES.joinpath(
            "regions", "gaps", "annotated",
            "{sample}",
            "{sample}.asm-{asm_unit}.{refgenome}.hprc-gaps.details.tsv"
        ),
        summary = DIR_RES.joinpath(
            "regions", "gaps", "annotated",
            "{sample}",
            "{sample}.asm-{asm_unit}.{refgenome}.hprc-gaps.summary.tsv"
        )
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt
    params:
        script=DIR_SCRIPTS.joinpath("eval_gaps", "add_qc_info.py"),
        asm_unit=lambda wildcards: f"asm-{wildcards.asm_unit}"
    shell:
        "{params.script} -g {input.gaps} -b {input.isect_table} "
        "-f {input.flagger} -n {input.nucfreq} "
        "-k {input.karyotypes} -c {input.chrom_assign} "
        "-s {wildcards.sample} -a {params.asm_unit} "
        "--out-table {output.table} --out-summary {output.summary}"


ERROR_THRESHOLD_MAP = {
    "1pct": 1.,
    "01pct": 0.1,
    "5pct": 5.
}


rule merge_hprc_gap_details:
    input:
        tables = expand(
            rules.annotate_gaps_with_qc_info.output.table,
            sample=SAMPLES,
            asm_unit=["hap1", "hap2"],
            refgenome=["t2tv2"]
        )
    output:
        summary = DIR_RES.joinpath(
            "regions", "gaps", "summary",
            "SAMPLES.vrk-ps-sseq.{refgenome}.hprc-gaps.ctg-summary-{err_t}.tsv"
        )
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    params:
        threshold = lambda wildcards: ERROR_THRESHOLD_MAP[wildcards.err_t],
    run:
        import pandas as pd
        import pathlib as pl

        merged = []
        sex_lut = dict()
        for table in sorted(input.tables):
            filename = pl.Path(table).name
            sample = filename.split(".")[0]
            au = filename.split(".")[2].split("-")[-1]

            df = pd.read_csv(table, sep="\t", header=0)
            df["aln_span"] = df["end"] - df["start"]
            hap_sex = df["sex"].iloc[0]
            sex_lut[f"{sample}.{au}"] = hap_sex

            select_flagger = df["flagger_pct"] < params.threshold
            select_nucfreq = df["nucfreq_pct"] < params.threshold
            select_covered = df["align_status"] == "covered"

            column_label = f"closed_at_err_{wildcards.err_t}"
            df[column_label] = 0
            select_rows = select_flagger & select_nucfreq & select_covered
            if select_rows.any():
                df.loc[select_rows, column_label] = 1

            select_open = df[column_label] == 0
            # same as: ~select_rows ...
            if select_open.any():
                df.loc[select_open, "seq"] = "no-cov"

            # for partial contig alignments, several rows may exist
            # w/ the same gap but, logically, this cannot represent
            # a fully covered gap, so drop those rows for deduplication purposes
            dups = df["gap_id"].duplicated(False)
            not_covered = df["align_status"] != "covered"
            # drop if partial and duplicated
            deselect = ~(dups & not_covered)
            df = df.loc[deselect, :].copy()

            dups = df["gap_id"].duplicated(False)
            if dups.any():
                # if still duplicates, this means that minimap and mashmap
                # generated two distinct alignments for the same sequence,
                # example: gap 09ab63f42e8d4d3810e203a3ea5a9cd2 in HG02106/H2
                # in this case, we select the alignment with largest span
                drop_indices = []
                for gap_id, alns in df.loc[dups, :].groupby("gap_id"):
                    max_span = alns["aln_span"].idxmax()
                    drop_indices.extend([i for i in alns.index if i != max_span])
                df.drop(drop_indices, axis=0, inplace=True)

            df = df[["chrom", "gap_start", "gap_end", "gap_id", "gap_length", "seq", column_label]]
            df.rename(
                {
                    "gap_start": "start",
                    "gap_end": "end",
                    "seq": f"{sample}.{au}.ctg",
                    column_label: f"{sample}.{au}.closed_{wildcards.err_t}"
                }, axis=1, inplace=True
            )
            df[f"{sample}.{au}.sex"] = hap_sex

            assert df["gap_id"].nunique() == df.shape[0]
            df.set_index(["chrom", "start", "end", "gap_id", "gap_length"], inplace=True)
            merged.append(df)

        merged = pd.concat(merged, axis=1, ignore_index=False, join="outer")
        ctg_columns = [c for c in merged.columns if c.endswith("ctg")]
        merged[ctg_columns] = merged[ctg_columns].fillna("no-ctg", inplace=False)

        stat_columns = [c for c in merged.columns if "closed" in c]
        merged[stat_columns] = merged[stat_columns].fillna(-1, inplace=False).astype(int)

        for hap_assm, hap_sex in sex_lut.items():
            merged[f"{hap_assm}.sex"] = merged[f"{hap_assm}.sex"].fillna(hap_sex, inplace=False)

        merged.sort_index(inplace=True)
        merged.to_csv(output.summary, sep="\t", header=True, index=True)
    # END OF RUN BLOCK


localrules: merge_hprc_gap_summaries
rule merge_hprc_gap_summaries:
    input:
        summaries = expand(
            rules.annotate_gaps_with_qc_info.output.summary,
            sample=SAMPLES,
            asm_unit=["hap1", "hap2"],
            allow_missing=True
        )
    output:
        summary = DIR_RES.joinpath(
            "regions", "gaps", "summary",
            "SAMPLES.vrk-ps-sseq.{refgenome}.hprc-gaps.summary.tsv"
        )
    run:
        import pandas as pd

        merged = pd.concat([
            pd.read_csv(summary_file, sep="\t", header=0)
            for summary_file in input.summaries
        ], axis=0, ignore_index=False)

        merged.sort_values(["sample", "asm_unit", "chrom", "max_stringency"], inplace=True)
        merged.to_csv(output.summary, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule run_all_annotate_gaps:
    input:
        tables = expand(
            rules.simplify_segdup_gap_intersection.output.tsv,
            sample=SAMPLES,
            asm_unit=MAIN_ASSEMBLY_UNITS,
            refgenome=["hg38", "t2tv2"],
            pct_id=["095", "098"]
        ),
        gaps = expand(
            rules.simplify_hprc_gap_intersection.output.tsv,
            sample=SAMPLES,
            asm_unit=MAIN_ASSEMBLY_UNITS,
            refgenome=["t2tv2"]
        ),
        norm_nf = expand(
            rules.filter_nucfreq_regions_to_bed.output,
            sample=SAMPLES
        ),
        split_flag = expand(
            rules.split_normalized_flagger.output,
            sample=SAMPLES
        ),
        gap_summaries = expand(
            rules.annotate_gaps_with_qc_info.output.summary,
            sample=SAMPLES,
            asm_unit=["hap1", "hap2"],
            refgenome=["t2tv2"]
        ),
        mrg_summary = expand(
            rules.merge_hprc_gap_summaries.output.summary,
            refgenome=["t2tv2"]
        ),
        ctg_summary = expand(
            rules.merge_hprc_gap_details.output.summary,
            refgenome=["t2tv2"],
            err_t=["1pct", "5pct", "01pct"]
        )
