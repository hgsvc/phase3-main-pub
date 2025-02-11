
rule filter_mosdepth_regions:
    input:
        bed = WORKDIR_EVAL.joinpath(
            "proc/50-postprocess/asm_ctg_readcov/mosdepth/{sample}.{read_type}.onlyPRI.{mapq}.wd",
            "{sample}.{read_type}.onlyPRI.{mapq}.regions.bed.gz"
        )
    output:
        bed_like = DIR_RES.joinpath(
            "asm_label_qc", "norm_tables",
            "read_depth", "{sample}.{read_type}.{mapq}.mosdepth-windowed.tsv.gz"
        )
    resources:
        mem_mb=lambda wildcards, attempt: 4096 * attempt
    run:
        import pandas as pd

        df = pd.read_csv(input.bed, sep="\t", header=None,
            names=["seq", "start", "end", "cov"]
        )
        df[["seqname", "tag"]] = df["seq"].str.rsplit(".", 1, expand=True)
        df = df.loc[df["tag"].isin(["hap1", "hap2", "unassigned"]), :].copy()
        df[["seqname", "start", "end", "cov"]].to_csv(
            output.bed_like, sep="\t", header=True, index=False
        )
    # END OF RUN BLOCK


rule run_all_filter_mosdepth_regions:
    input:
        tables = expand(
            rules.filter_mosdepth_regions.output.bed_like,
            sample=SAMPLES,
            read_type=["hifi", "ont"],
            mapq=["mq00", "mq60"]
        )
