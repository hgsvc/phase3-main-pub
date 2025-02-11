
localrules: split_cnv_table_by_sample
rule split_cnv_table_by_sample:
    input:
        bed = ILLCNV_ANNOTATION
    output:
        bed = expand(
            DIR_PROC.joinpath(
                "regions", "roi", "illcnv",
                "{sample}.illcnv-calls.hg38.bed"
            ),
            sample=SAMPLES
        )
    run:
        import pandas as pd

        assert ASSEMBLER == "verkko"

        cnvs = pd.read_csv(input.bed, sep="\t", header=0)
        done_samples = set()
        out_header = None
        for sample, subset in cnvs.groupby("sample"):
            out_path = DIR_PROC.joinpath(
                "regions", "roi", "illcnv",
                f"{sample}.vrk-ps-sseq.illcnv-calls.hg38.bed"
            )
            subset.to_csv(out_path, sep="\t", header=True, index=False)
            out_header = subset.columns.values
            done_samples.add(sample)

        remainder = set(PLAIN_SAMPLES) - done_samples
        out_header = "\t".join(out_header)

        for sample in remainder:
            out_path = DIR_PROC.joinpath(
                "regions", "roi", "illcnv",
                f"{sample}.vrk-ps-sseq.illcnv-calls.hg38.bed"
            )
            with open(out_path, "w") as mock:
                _ = mock.write(out_header + "\n")

            with open(out_path.with_suffix(".bed.MOCK"), "w"):
                pass

    # END OF RUN BLOCK


rule intersect_cnvs_with_alnblocks:
    input:
        cnv_bed = DIR_PROC.joinpath(
            "regions", "roi", "illcnv",
            "{sample}.illcnv-calls.hg38.bed"
        ),
        hap1_blocks = WORKDIR_EVAL.joinpath(
            "results", "regions", "{sample}",
            "{sample}.asm-hap1.hg38.ctg-aln-gap.ref-coord.bed"
        ),
        hap2_blocks = WORKDIR_EVAL.joinpath(
            "results", "regions", "{sample}",
            "{sample}.asm-hap2.hg38.ctg-aln-gap.ref-coord.bed"
        ),
        un_blocks = WORKDIR_EVAL.joinpath(
            "results", "regions", "{sample}",
            "{sample}.asm-unassigned.hg38.ctg-aln-gap.ref-coord.bed"
        )
    output:
        isect_hap1 = DIR_PROC.joinpath(
            "regions", "roi", "illcnv", "intersect",
            "{sample}.asm-hap1.hg38.ctg-aln-gap.illcnv.tsv"
        ),
        isect_hap2 = DIR_PROC.joinpath(
            "regions", "roi", "illcnv", "intersect",
            "{sample}.asm-hap2.hg38.ctg-aln-gap.illcnv.tsv"
        ),
        isect_un = DIR_PROC.joinpath(
            "regions", "roi", "illcnv", "intersect",
            "{sample}.asm-unassigned.hg38.ctg-aln-gap.illcnv.tsv"
        ),
    conda:
        DIR_ENVS.joinpath("bedtools.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt
    shell:
        "bedtools intersect -wao -a {input.cnv_bed} -b {input.hap1_blocks} > {output.isect_hap1}"
            " && "
        "bedtools intersect -wao -a {input.cnv_bed} -b {input.hap2_blocks} > {output.isect_hap2}"
            " && "
        "bedtools intersect -wao -a {input.cnv_bed} -b {input.un_blocks} > {output.isect_un}"


rule run_all_illcnv:
    input:
        isects = expand(
            DIR_PROC.joinpath(
                "regions", "roi", "illcnv", "intersect",
                "{sample}.asm-{asm_unit}.hg38.ctg-aln-gap.illcnv.tsv"
            ),
            sample=SAMPLES,
            asm_unit=MAIN_ASSEMBLY_UNITS
        )
