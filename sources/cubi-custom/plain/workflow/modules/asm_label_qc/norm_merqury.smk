
rule merge_merqury_kmer_tracks:
    input:
        bed_files = collect_merqury_kmer_tracks
    output:
        bed_like = temp(DIR_PROC.joinpath(
            "asm_label_qc", "norm_merqury",
            "{sample}.asmonly-kmer.bed"
        ))
    conda:
        DIR_ENVS.joinpath("bedtools.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    shell:
        "cat {input.bed_files} | sort -V -k1,1 -k2,2n"
            " | "
        "bedtools merge -i /dev/stdin > {output.bed_like}"


rule normalize_merqury_kmer_tracks:
    input:
        bed_like = rules.merge_merqury_kmer_tracks.output.bed_like
    output:
        bed_like = DIR_RES.joinpath(
            "asm_label_qc", "norm_tables",
            "merqury", "{sample}.merqury-asmonly-kmer.tsv.gz"
        )
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    run:
        import pandas as pd
        df = pd.read_csv(input.bed_like, sep="\t", header=None, names=["chrom", "start", "end"])
        df["name"] = "merqury_asmonly_kmer"
        df.to_csv(output.bed_like, sep="\t", header=True, index=False)


localrules: add_merqury_merge_label
rule add_merqury_merge_label:
    input:
        bed = rules.normalize_merqury_kmer_tracks.output.bed_like
    output:
        bed = DIR_RES.joinpath(
            "asm_label_qc", "merge_tables", "by-sample",
            "{sample}", "{sample}.merqury.mrg-labels.bed"
        )
    run:
        import pandas as pd
        df = pd.read_csv(input.bed, sep="\t", header=0)
        df["raw_label"] = "MRQURY"
        df["length"] = (df["end"] - df["start"]).astype(int)
        df["merge_label"] = df["raw_label"] + "::" + df["length"].astype(str)
        df = df[["chrom", "start", "end", "merge_label"]]
        df.to_csv(output.bed, sep="\t", header=False, index=False)
    # END OF RUN BLOCK


rule run_normalize_merqury_results:
    input:
        tables = expand(
            rules.normalize_merqury_kmer_tracks.output.bed_like,
            sample=SAMPLES
        ),
        mrg_bed = expand(
            rules.add_merqury_merge_label.output.bed,
            sample=SAMPLES
        )
