
localrules: normalize_merqury_qv_estimates
rule normalize_merqury_qv_estimates:
    """Change behavior: compute global / whole-genome
    estimates as weighted averages over per-sequence estimates
    in rule downstream
    """
    input:
        detail_file = lambda wildcards: find_merqury_output_file(
            wildcards.sample, "qv_detail", phased_only=(wildcards.assembler != "verkko")
        )
    output:
        tsv = DIR_PROC.joinpath(
            "merqury", "{assembler}", "{sample}",
            "{sample}.merqury-qv-est.per-seq.tsv"
        )
    run:
        import pandas as pd
        import numpy as np

        def assign_asm_unit(seqname):

            if any(h in seqname for h in ["h1", "hap1", "haplotype1"]):
                return "hap1"
            elif any(h in seqname for h in ["h2", "hap2", "haplotype2"]):
                return "hap2"
            elif "unassigned" in seqname:
                return "unassigned"
            else:
                raise ValueError(seqname)

        concat = []
        for hap_file in input.detail_file:
            detail = pd.read_csv(
                hap_file, sep="\t", header=None,
                names=["sequence", "error_bp", "total_adj_bp", "qv_est", "error_rate"]
            )
            detail["sample"] = wildcards.sample
            detail["asm_unit"] = detail["sequence"].apply(assign_asm_unit)
            no_error = np.isinf(detail["qv_est"].values)
            detail.loc[no_error, "qv_est"] = 99
            concat.append(detail)
        concat = pd.concat(concat, axis=0, ignore_index=False)

        concat = concat[["sample", "asm_unit", "sequence", "error_bp", "total_adj_bp", "qv_est", "error_rate"]]
        with open(output.tsv, "w") as table:
            concat.to_csv(table, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


localrules: compute_global_merqury_qv_estimates
rule compute_global_merqury_qv_estimates:
    input:
        table = rules.normalize_merqury_qv_estimates.output.tsv
    output:
        tsv = DIR_RES.joinpath(
            "merqury", "{assembler}", "{sample}",
            "{sample}.merqury-qv-est.tsv"
        )
    run:
        import pandas as pd
        import numpy as np

        df = pd.read_csv(input.table, sep="\t", header=0)

        combinations = [
            ("hap1",), ("hap2",), ("hap1", "hap2"),
            ("hap1", "hap2", "unassigned"), ("unassigned",)
        ]
        labels = ["hap1", "hap2", "phased", "wg", "unassigned"]

        summary_df = []
        for comb, label in zip(combinations, labels):
            subset = df.loc[df["asm_unit"].isin(comb), :]
            if subset.empty:
                continue
            total_error = subset["error_bp"].sum()
            total_adj_len = subset["total_adj_bp"].sum()
            qv_est = np.average(subset["qv_est"].values, weights=subset["total_adj_bp"].values)
            error_rate = np.average(subset["error_rate"].values, weights=subset["total_adj_bp"].values)
            summary_df.append(
                [subset["sample"].iloc[0], label, "total", total_error, total_adj_len, qv_est, error_rate]
            )
        summary_df = pd.DataFrame(
            summary_df,
            columns=["sample", "asm_unit", "sequence", "error_bp", "total_adj_bp", "qv_est", "error_rate"]
        )

        df = pd.concat([df, summary_df], axis=0, ignore_index=False)
        df.sort_values(["asm_unit", "total_adj_bp"], ascending=[True, False], inplace=True)

        df.to_csv(output.tsv, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule merge_merqury_normalized_qv_estimates:
    input:
        tables = expand(
            rules.compute_global_merqury_qv_estimates.output.tsv,
            sample=SAMPLES,
            allow_missing=True
        )
    output:
        tsv = DIR_RES.joinpath(
            "merqury", "{assembler}.qv-est.tsv"
        )
    run:
        import pandas as pd
        import io

        merge = []
        for table_file in input.tables:
            table_buffer = io.StringIO()
            with open(table_file, "r") as table:
                for line in table:
                    if line.startswith("#"):
                        continue
                    table_buffer.write(line)
            table_buffer.seek(0)
            df = pd.read_csv(table_buffer, sep="\t", header=0)
            merge.append(df)
        merge = pd.concat(merge, axis=0, ignore_index=False)
        merge.sort_values(["sample", "asm_unit", "sequence"], inplace=True)
        merge.to_csv(output.tsv, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


################################
### ABOVE: QV ESTIMATES
################################
### BELOW: K-MER COMPLETENESS
################################

localrules: normalize_merqury_kmer_completeness_hifiasm
rule normalize_merqury_kmer_completeness_hifiasm:
    input:
        table = lambda wildcards: find_merqury_output_file(wildcards.sample, "completeness")
    output:
        tsv = DIR_RES.joinpath(
            "merqury", "{assembler}", "{sample}",
            "{sample}.merqury-kmer-completeness.tsv"
        )
    wildcard_constraints:
        assembler="hifiasm"
    run:
        import pandas as pd

        rows = []
        with open(input.table[0], "r") as listing:
            for line in listing:
                columns = line.strip().split()
                assert columns[1] == "all"  # no clue ...
                pct_complete = float(columns[-1])
                kmer_total = int(columns[-2])
                kmer_found = int(columns[-3])
                asm_unit = None
                if "hap1" in columns[0]:
                    asm_unit = "hap1"
                elif "hap2" in columns[0]:
                    asm_unit = "hap2"
                elif "both" in columns[0].lower():
                    asm_unit = "wg"
                else:
                    raise ValueError(line.strip())
                assert asm_unit is not None
                rows.append(
                    (wildcards.sample, asm_unit, kmer_found, kmer_total, pct_complete)
                )
                if asm_unit == "wg":
                    rows.append(
                        (wildcards.sample, "phased", kmer_found, kmer_total, pct_complete)
                    )
        df = pd.DataFrame.from_records(
            rows, columns=["sample", "asm_unit", "kmer_present", "kmer_total", "kmer_present_pct"]
        )
        with open(output.tsv, "w") as dump:
            df.to_csv(dump, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


localrules: normalize_merqury_kmer_completeness_verkko
rule normalize_merqury_kmer_completeness_verkko:
    input:
        phased = lambda wildcards: find_merqury_output_file(wildcards.sample, "completeness", phased_only=True),
        all_incl = lambda wildcards: find_merqury_output_file(wildcards.sample, "completeness", phased_only=False)
    output:
        tsv = DIR_RES.joinpath(
            "merqury", "{assembler}", "{sample}",
            "{sample}.merqury-kmer-completeness.tsv"
        )
    wildcard_constraints:
        assembler="verkko"
    run:
        import pandas as pd

        with open(input.phased[0], "r") as listing:
            asm_hap1, _, hap1_present, hap1_total, hap1_pct = listing.readline().strip().split()
            asm_hap2, _, hap2_present, hap2_total, hap2_pct = listing.readline().strip().split()
            asm_phased, _, phased_present, phased_total, phased_pct = listing.readline().strip().split()

        with open(input.all_incl[0], "r") as listing:
            asm_incl, _, incl_present, incl_total, incl_pct = listing.readline().strip().split()
            assert "unassigned" in asm_incl.lower()
            asm_hap2, _, hap2_present, hap2_total, hap2_pct = listing.readline().strip().split()
            wg, _, wg_present, wg_total, wg_pct = listing.readline().strip().split()

        # DEBUG: waiting for William H. to shed light on this / mismatch for HG00733
        if wildcards.sample != "HG00733.vrk-ps-sseq":
            assert int(hap1_total) == int(incl_total)  # check same k-mer db
        else:
            # resolve / confirm this
            hap1_total = int(incl_total)
            hap2_total = int(incl_total)
            phased_total = int(incl_total)
        unassigned_present = int(incl_present) - int(hap1_present)
        unassigned_total = int(incl_total)
        unassigned_pct = float(unassigned_present / unassigned_total)

        columns=["sample", "asm_unit", "kmer_present", "kmer_total", "kmer_present_pct"]

        rows = [
            [wildcards.sample, "wg", int(wg_present), int(wg_total), float(wg_pct)],
            [wildcards.sample, "phased", int(phased_present), int(phased_total), float(phased_pct)],
            [wildcards.sample, "hap1", int(hap1_present), int(hap1_total), float(hap1_pct)],
            [wildcards.sample, "hap2", int(hap2_present), int(hap2_total), float(hap2_pct)],
            [wildcards.sample, "unassigned", int(unassigned_present), int(unassigned_total), round(float(unassigned_pct), 4)]
        ]

        df = pd.DataFrame.from_records(rows, columns=columns)

        with open(output.tsv, "w") as dump:
            #_ = dump.write(f"# {shorten_merqury_file_path(input.table[0])}\n")
            df.to_csv(dump, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule merge_merqury_normalized_kmer_completeness:
    input:
        tables = expand(
            DIR_RES.joinpath(
                "merqury", "{assembler}", "{sample}",
                "{sample}.merqury-kmer-completeness.tsv"
            ), sample=SAMPLES, allow_missing=True
        )
    output:
        tsv = DIR_RES.joinpath(
            "merqury", "{assembler}.kmer-completeness.tsv"
        )
    run:
        import pandas as pd

        merge = []
        for table in input.tables:
            df = pd.read_csv(table, sep="\t", header=0, comment="#")
            merge.append(df)
        merge = pd.concat(merge, axis=0, ignore_index=False)
        merge.sort_values(["sample", "asm_unit"], inplace=True)
        merge.to_csv(output.tsv, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


################################
### ABOVE: K-MER COMPLETENESS
################################
### BELOW: K-MER BED TRACKS
################################


rule merge_merqury_asm_only_kmers_hifiasm:
    input:
        bed_files = lambda wildcards: find_merqury_output_file(wildcards.sample, "errors")
    output:
        bed_hap1 = DIR_RES.joinpath(
            "merqury", "{assembler}", "{sample}",
            "{sample}.asm-hap1.asm-only-kmer.bed.gz"
        ),
        bed_hap2 = DIR_RES.joinpath(
            "merqury", "{assembler}", "{sample}",
            "{sample}.asm-hap2.asm-only-kmer.bed.gz"
        ),
        tbi_hap1 = DIR_RES.joinpath(
            "merqury", "{assembler}", "{sample}",
            "{sample}.asm-hap1.asm-only-kmer.bed.gz.tbi"
        ),
        tbi_hap2 = DIR_RES.joinpath(
            "merqury", "{assembler}", "{sample}",
            "{sample}.asm-hap2.asm-only-kmer.bed.gz.tbi"
        )
    wildcard_constraints:
        assembler="hifiasm"
    conda:
        DIR_ENVS.joinpath("bedtools.yaml")
    params:
        bed_hap1 = lambda wildcards, input: input.bed_files[0],
        bed_hap2 = lambda wildcards, input: input.bed_files[1]
    shell:
        "bedtools merge -i {params.bed_hap1} | bgzip > {output.bed_hap1} && tabix -p bed {output.bed_hap1}"
            " && "
        "bedtools merge -i {params.bed_hap2} | bgzip > {output.bed_hap2} && tabix -p bed {output.bed_hap2}"


rule merge_merqury_asm_only_kmers_verkko:
    input:
        bed_files = lambda wildcards: find_merqury_output_file(wildcards.sample, "errors", phased_only=False)
    output:
        tmp_hap1 = temp(
            DIR_PROC.joinpath("merqury", "{assembler}", "{sample}.asm-hap1.kmer-tmp.bed")
        ),
        tmp_unassign = temp(
            DIR_PROC.joinpath("merqury", "{assembler}", "{sample}.asm-unassign.kmer-tmp.bed")
        ),
        bed_hap1 = DIR_RES.joinpath(
            "merqury", "{assembler}", "{sample}",
            "{sample}.asm-hap1.asm-only-kmer.bed.gz"
        ),
        bed_hap2 = DIR_RES.joinpath(
            "merqury", "{assembler}", "{sample}",
            "{sample}.asm-hap2.asm-only-kmer.bed.gz"
        ),
        bed_unassign = DIR_RES.joinpath(
            "merqury", "{assembler}", "{sample}",
            "{sample}.asm-unassigned.asm-only-kmer.bed.gz"
        ),
        tbi_hap1 = DIR_RES.joinpath(
            "merqury", "{assembler}", "{sample}",
            "{sample}.asm-hap1.asm-only-kmer.bed.gz.tbi"
        ),
        tbi_hap2 = DIR_RES.joinpath(
            "merqury", "{assembler}", "{sample}",
            "{sample}.asm-hap2.asm-only-kmer.bed.gz.tbi"
        ),
        tbi_unassign = DIR_RES.joinpath(
            "merqury", "{assembler}", "{sample}",
            "{sample}.asm-unassigned.asm-only-kmer.bed.gz.tbi"
        )
    wildcard_constraints:
        assembler="verkko"
    conda:
        DIR_ENVS.joinpath("bedtools.yaml")
    params:
        bed_hap1 = lambda wildcards, input: input.bed_files[0],
        bed_hap2 = lambda wildcards, input: input.bed_files[1]
    shell:
        "grep -F haplotype1 {params.bed_hap1} > {output.tmp_hap1}"
            " && "
        "bedtools merge -i {output.tmp_hap1} | bgzip > {output.bed_hap1} && tabix -p bed {output.bed_hap1}"
            " && "
        "grep -F unassigned {params.bed_hap1} > {output.tmp_unassign}"
            " && "
        "bedtools merge -i {output.tmp_unassign} | bgzip > {output.bed_unassign} && tabix -p bed {output.bed_unassign}"
            " && "
        "bedtools merge -i {params.bed_hap2} | bgzip > {output.bed_hap2} && tabix -p bed {output.bed_hap2}"


rule run_all_merqury_postprocess:
    input:
        qv_est = expand(
            rules.merge_merqury_normalized_qv_estimates.output.tsv,
            assembler=[ASSEMBLER]
        ),
        kmer_completeness = expand(
            rules.merge_merqury_normalized_kmer_completeness.output.tsv,
            assembler=[ASSEMBLER]
        ),
        bed_files = expand(
            DIR_RES.joinpath(
                "merqury", "{assembler}", "{sample}",
                "{sample}.asm-{asm_unit}.asm-only-kmer.bed.gz"
            ),
            sample=SAMPLES,
            assembler=[ASSEMBLER],
            asm_unit=MAIN_ASSEMBLY_UNITS
        )
