
rule bin_nucfreq_regions_by_coverage:
    input:
        bed = WORKDIR_EVAL.joinpath(
            "results/regions/{sample}",
            "{sample}.nucfreq.covann.tsv.gz"
        )
    output:
        bed = DIR_RES.joinpath(
            "asm_label_qc", "norm_tables", "nucfreq",
            "{sample}.nucfreq-cov-bin.tsv.gz"
        )
    run:
        import pandas as pd
        import sys

        df = pd.read_csv(input.bed, sep="\t", header=0)
        bins = [0, 50, 90, 110, 150, 200, sys.maxsize]
        df["score"] = pd.cut(
            df["hifi_pct_median_cov"].values,
            bins, right=False, include_lowest=True,
            labels=False
        )
        df = df[["contig", "start", "end", "score", "asm_unit", "num_hets"]]
        df.to_csv(output.bed, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


localrules: binarize_nucfreq_output
rule binarize_nucfreq_output:
    input:
        bed = rules.bin_nucfreq_regions_by_coverage.output.bed
    output:
        bed = DIR_RES.joinpath(
            "asm_label_qc", "norm_tables", "nucfreq",
            "{sample}.nucfreq-binary.tsv.gz"
        )
    run:
        import pandas as pd

        df = pd.read_csv(input.bed, sep="\t", header=0, usecols=["contig", "start", "end"])
        df.to_csv(output.bed, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


localrules: add_nucfreq_merge_label
rule add_nucfreq_merge_label:
    input:
        bed = rules.binarize_nucfreq_output.output.bed
    output:
        bed = DIR_RES.joinpath(
            "asm_label_qc", "merge_tables", "by-sample",
            "{sample}", "{sample}.nucfreq.mrg-labels.bed"
        )
    run:
        import pandas as pd
        df = pd.read_csv(input.bed, sep="\t", header=0)
        df["raw_label"] = "NUCFRQ"
        df["length"] = (df["end"] - df["start"]).astype(int)
        df["merge_label"] = df["raw_label"] + "::" + df["length"].astype(str)
        df = df[["contig", "start", "end", "merge_label"]]
        df.to_csv(output.bed, sep="\t", header=False, index=False)
    # END OF RUN BLOCK


rule run_all_bin_nucfreq_regions:
    input:
        beds = expand(
            rules.bin_nucfreq_regions_by_coverage.output.bed,
            sample=SAMPLES
        ),
        mrg_bed = expand(
            rules.add_nucfreq_merge_label.output.bed,
            sample=SAMPLES
        )
