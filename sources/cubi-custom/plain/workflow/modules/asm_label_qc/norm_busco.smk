
localrules: normalize_busco_issue_annotation
rule normalize_busco_issue_annotation:
    input:
        bed = WORKDIR_EVAL.joinpath(
            "results/regions", "{sample}",
            "{sample}.busco.primates_odb10.issues.bed"
        )
    output:
        bed = DIR_RES.joinpath(
            "asm_label_qc", "norm_tables", "busco",
            "{sample}.busco-issues.bed"
        )
    run:
        import pandas as pd
        df = pd.read_csv(input.bed, sep="\t", header=0)
        df.drop("gene", axis=1, inplace=True)
        relabel = {
            "Fragmented": 1,
            "Duplicated": 2
        }
        df.insert(4, "score", 1)
        df["score"] = df["label"].replace(relabel, inplace=False)
        df["score"] = df["score"].astype(int)  # capture other potential labels
        df.to_csv(output.bed, sep="\t", header=True, index=False)
    # END OF RUN BLOCK



localrules: add_busco_merge_label
rule add_busco_merge_label:
    input:
        bed = rules.normalize_busco_issue_annotation.output.bed
    output:
        bed = DIR_RES.joinpath(
            "asm_label_qc", "merge_tables", "by-sample",
            "{sample}", "{sample}.busco.mrg-labels.bed"
        )
    run:
        import pandas as pd
        df = pd.read_csv(input.bed, sep="\t", header=0)
        rename = {
            "Duplicated": "BSCDUP",
            "Fragmented": "BSCFRG"
        }
        df["raw_label"] = df["label"].replace(rename, inplace=False)
        df["length"] = (df["end"] - df["start"]).astype(int)
        df["merge_label"] = df["raw_label"] + "::" + df["length"].astype(str)

        df = df[["#seq_name", "start", "end", "merge_label"]]
        df.to_csv(output.bed, sep="\t", header=False, index=False)
    # END OF RUN BLOCK


rule run_all_normalize_busco:
    input:
        bed = expand(
            rules.normalize_busco_issue_annotation.output.bed,
            sample=SAMPLES
        ),
        mrg = expand(
            rules.add_busco_merge_label.output.bed,
            sample=SAMPLES
        )

