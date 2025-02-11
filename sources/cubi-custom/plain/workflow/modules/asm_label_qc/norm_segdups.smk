
rule split_segdup_annotation:
    input:
        folder = SEGDUP_ROOT_FOLDER
    output:
        sd95 = DIR_RES.joinpath(
            "asm_label_qc", "norm_tables",
            "segdups", "{sample}.sd-095.tsv.gz"
        ),
        sd98 = DIR_RES.joinpath(
            "asm_label_qc", "norm_tables",
            "segdups", "{sample}.sd-098.tsv.gz"
        ),
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    run:
        import pandas as pd

        sample_id, assm_id = wildcards.sample.split(".", 1)
        if assm_id == "vrk-ps-sseq":
            assembler = "verkko"
        elif assm_id == "hsm-ps-sseq":
            assembler = "hifiasm"
        else:
            raise ValueError(assm_id)

        bed_files = sorted(SEGDUP_ROOT_FOLDER.joinpath(
            f"{assembler}"
        ).glob(f"{sample_id}*.bed"))
        assert len(bed_files) == 2

        sd95 = []
        sd98 = []
        for bed_file in bed_files:
            df = pd.read_csv(
                bed_file, sep="\t", header=0,
                usecols=["#chr1", "start1", "end1", "fracMatch"]
            )
            df.rename({"#chr1": "seq", "start1": "start", "end1": "end"}, axis=1, inplace=True)

            lower = df["fracMatch"] > 0.95
            upper = df["fracMatch"] <= 0.98
            sub1 = df.loc[lower & upper, ["seq", "start", "end"]].copy()
            sd95.append(sub1)

            lower = df["fracMatch"] > 0.98
            upper = df["fracMatch"] <= 1.
            sub2 = df.loc[lower & upper, ["seq", "start", "end"]].copy()
            sd98.append(sub2)

        sd95 = pd.concat(sd95, axis=0, ignore_index=False)
        sd95.sort_values(["seq", "start"], inplace=True)
        sd95.drop_duplicates(inplace=True)
        sd95.to_csv(output.sd95, header=True, index=False, sep="\t")

        sd98 = pd.concat(sd98, axis=0, ignore_index=False)
        sd98.sort_values(["seq", "start"], inplace=True)
        sd98.drop_duplicates(inplace=True)
        sd98.to_csv(output.sd98, header=True, index=False, sep="\t")

    # END OF RUN BLOCK


localrules: add_segdup_merge_label
rule add_segdup_merge_label:
    input:
        sd95 = DIR_RES.joinpath(
            "asm_label_qc", "norm_tables",
            "segdups", "{sample}.sd-095.tsv.gz"
        ),
        sd98 = DIR_RES.joinpath(
            "asm_label_qc", "norm_tables",
            "segdups", "{sample}.sd-098.tsv.gz"
        ),
    output:
        bed = DIR_RES.joinpath(
            "asm_label_qc", "merge_tables", "by-sample",
            "{sample}", "{sample}.segdup.mrg-labels.bed"
        )
    run:
        import pandas as pd

        concat = []
        for input_file in [input.sd95, input.sd98]:
            if "sd-095" in str(input_file):
                label = "SEGD95"
            elif "sd--098" in str(input_file):
                label = "SEGD98"
            else:
                raise

            df = pd.read_csv(input_file, sep="\t", header=0)
            df["raw_label"] = label
            df["length"] = (df["end"] - df["start"]).astype(int)
            df["merge_label"] = df["raw_label"] + "::" + df["length"].astype(str)
            df = df[["chrom", "start", "end", "merge_label"]]
            concat.append(df)
        concat = pd.concat(concat, axis=0, ignore_index=False)
        concat.sort_values(["seq", "start", "end"], inplace=True)
        concat.to_csv(output.bed, sep="\t", header=False, index=False)
    # END OF RUN BLOCK


rule run_all_normalize_segdups:
    input:
        segdups = expand(
            rules.split_segdup_annotation.output,
            sample=SAMPLES
        ),
        mrg_bed = expand(
            rules.add_segdup_merge_label.output.bed,
            sample=SAMPLES
        )
