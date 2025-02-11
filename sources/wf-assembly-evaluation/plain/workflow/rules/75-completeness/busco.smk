
localrules: compleasm_list_available_busco_lineages
rule compleasm_list_available_busco_lineages:
    """NB: this makes use of the wildcard to allow for later
    db updates w/o touching all of the already generated outputs
    """
    input:
        busco_db = DIR_GLOBAL_REF.joinpath("busco_db")
    output:
        lineage = DIR_PROC.joinpath(
            "75-completeness", "busco", "lineage_{odb_name}.check.txt"
        )
    conda:
        DIR_ENVS.joinpath("biotools", "compleasm.yaml")
    shell:
        "compleasm list --local --library_path {input.busco_db} > {output.lineage}"


localrules: confirm_busco_lineage_exists
rule confirm_busco_lineage_exists:
    input:
        lineage = rules.compleasm_list_available_busco_lineages.output.lineage
    output:
        exists = DIR_PROC.joinpath(
            "75-completeness", "busco", "lineage_{odb_name}.exists.txt"
        )
    run:
        lineage_found = False
        with open(input.lineage, "r") as listing:
            first_line = listing.readline()
            if not first_line.startswith("Local available lineages"):
                err_msg = (
                    "Unexpected first line in compleasm / busco lineage file:\n"
                    f"{input.lineage}\n"
                    f"First line: {first_line.strip()}"
                )
                logerr(err_msg)
                raise ValueError(err_msg)
            for line in listing:
                if wildcards.odb_name in line or line.strip() in wildcards.odb_name:
                    lineage_found = True
                    break
        if not lineage_found:
            err_msg = "Requested BUSCO lineage not locally available: {wildcards.odb_name}"
            logerr(err_msg)
            raise ValueError(err_msg)
        with open(output.exists, "w") as dump:
            _ = dump.write(get_timestamp() + "\n")
            _ = dump.write(wildcards.odb_name + "\n")
    # END OF RUN BLOCK


rule compleasm_busco_mode:
    input:
        asm = get_asm_unit,
        lineage_exists = rules.confirm_busco_lineage_exists.output.exists,
        busco_db = DIR_GLOBAL_REF.joinpath("busco_db")
    output:
        summary = DIR_PROC.joinpath(
            "75-completeness", "busco", "{sample}.{asm_unit}.{odb_name}.wd",
            "summary.txt"
        )
    log:
        DIR_LOG.joinpath("75-completeness", "busco", "{sample}.{asm_unit}.{odb_name}.compleasm.log")
    benchmark:
        DIR_RSRC.joinpath("75-completeness", "busco", "{sample}.{asm_unit}.{odb_name}.compleasm.log")
    conda:
        DIR_ENVS.joinpath("biotools", "compleasm.yaml")
    threads: CPU_LOW
    resources:
        mem_mb = lambda wildcards, attempt: 32768 * attempt,
        time_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        wd=lambda wildcards, output: pathlib.Path(output.summary).parent,
    shell:
        "compleasm run --mode busco -L {input.busco_db} -l {wildcards.odb_name} "
        "--threads {threads} -o {params.wd} -a {input.asm} &> {log}"


localrules: normalize_compleasm_summary
rule normalize_compleasm_summary:
    input:
        txt = rules.compleasm_busco_mode.output.summary
    output:
        tsv = DIR_PROC.joinpath(
            "75-completeness", "busco", "{sample}.{asm_unit}.{odb_name}.compleasm-summary.tsv",
        )
    wildcard_constraints:
        asm_unit="(" + "|".join(ASSEMBLY_UNITS_MAIN) + ")"
    run:
        stats_map = {
            "S": "singleton", "D": "duplicated", "F": "fragmented",
            "I": "interspersed", "M": "missing", "N": "total"
        }
        data_row = [wildcards.sample, wildcards.asm_unit, wildcards.odb_name]
        header_row = ["sample", "asm_unit", "odb_name"]
        with open(input.txt, "r") as text:
            for line in text:
                if line.startswith("#"):
                    assert wildcards.odb_name in line
                    continue
                if line.startswith("N"):
                    percentage = 100.
                    statistic, count = line.strip().split(":")
                else:
                    statistic, count = line.strip().split(",")
                    count = int(count.strip())
                    statistic, percentage = statistic.split(":")
                    percentage = float(percentage.strip("%"))
                stat_name = stats_map[statistic]
                header_row.append(stat_name + "_count")
                data_row.append(str(count))
                header_row.append(stat_name + "_pct")
                data_row.append(str(percentage))
        with open(output.tsv, "w") as table:
            _ = table.write("\t".join(header_row) + "\n")
            _ = table.write("\t".join(data_row) + "\n")
    # END OF RUN BLOCK


localrules: merge_compleasm_output_tables
rule merge_compleasm_output_tables:
    input:
        summaries = expand(
            rules.compleasm_busco_mode.output.summary,
            asm_unit=ASSEMBLY_UNITS_MAIN,
            allow_missing=True
        )
    output:
        tsv = DIR_RES.joinpath(
            "reports", "completeness", "busco_per_sample",
            "{sample}.busco.{odb_name}.merged.tsv.gz"
        ),
        issues = DIR_RES.joinpath(
            "regions", "{sample}",
            "{sample}.busco.{odb_name}.issues.bed"
        )
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt
    run:
        import pandas as pd
        import pathlib as pl

        table_header = [
            "gene", "label", "seq_name",
            "start", "end", "strand",
            "score", "length", "identity",
            "fraction", "frameshift_events",
            "best_gene", "codons"
        ]

        merged = []
        all_issues = []
        for summary_file in input.summaries:
            summ_file_path = pl.Path(summary_file)
            file_parts = summ_file_path.parent.name.rsplit(".", 3)
            assert file_parts[0] == wildcards.sample, file_parts
            asm_unit = file_parts[1]
            label_column = f"{asm_unit}_label"
            folder = pl.Path(summary_file).parent
            table_file = folder.joinpath(
                f"{wildcards.odb_name}", "full_table.tsv"
            )
            assert table_file.is_file()
            df = pd.read_csv(table_file, sep="\t", header=None, skiprows=1, names=table_header)
            to_merge = df[table_header[:2]].copy()
            to_merge.rename({"label": label_column}, axis=1, inplace=True)
            to_merge.drop_duplicates("gene", inplace=True)
            to_merge.set_index("gene", inplace=True)
            merged.append(to_merge)
            has_start = ~pd.isnull(df["start"])
            issues = df.loc[(df["label"] != "Single") & has_start, :].copy()
            issues["asm_unit"] = asm_unit
            issues = issues[["seq_name", "start", "end", "label", "gene", "asm_unit"]].copy()
            issues["start"] = issues["start"].astype(int)
            issues["end"] = issues["end"].astype(int)
            all_issues.append(issues)

        merged = pd.concat(merged, axis=1, ignore_index=False)
        merged.sort_index(inplace=True)
        merged.to_csv(output.tsv, sep="\t", header=True, index=True, index_label="gene")

        all_issues = pd.concat(all_issues, axis=0, ignore_index=False)
        all_issues.sort_values(["seq_name", "start", "end"], inplace=True)

        with open(output.issues, "w") as dump:
            _ = dump.write("#")
            all_issues.to_csv(dump, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


localrules: aggregate_compleasm_summaries
rule aggregate_compleasm_summaries:
    """URGENT TODO --- odb names as parameter!
    """
    input:
        tsv = expand(
            rules.normalize_compleasm_summary.output.tsv,
            sample=SAMPLES,
            asm_unit=ASSEMBLY_UNITS_MAIN,
            odb_name=["eukaryota_odb10", "primates_odb10"]
        )
    output:
        tsv = DIR_RES.joinpath(
            "reports", "completeness", "busco-summary{run_id}.tsv"
        )
    run:
        import pandas as pd

        concat = []
        for table_file in input.tsv:
            df = pd.read_csv(table_file, sep="\t", header=0)
            concat.append(df)
        concat = pd.concat(concat, axis=0, ignore_index=False)
        concat.sort_values(["sample", "asm_unit"], inplace=True)
        concat.to_csv(output.tsv, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


# TODO: make odb db name parameter
rule run_all_compleasm:
    input:
        report = expand(
            rules.aggregate_compleasm_summaries.output.tsv,
            run_id=RUN_SUFFIX
        ),
        sample_merge = expand(
            rules.merge_compleasm_output_tables.output.tsv,
            sample=SAMPLES,
            odb_name=["eukaryota_odb10", "primates_odb10"]
        ),
        issues = expand(
            rules.merge_compleasm_output_tables.output.issues,
            sample=SAMPLES,
            odb_name=["eukaryota_odb10", "primates_odb10"]
        )




