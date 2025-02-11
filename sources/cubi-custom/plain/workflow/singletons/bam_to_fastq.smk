import pandas as pd

MISSING_FASTQ_FILE = config["missing"]

MISSING_FASTQ = pd.read_csv(
    MISSING_FASTQ_FILE,
    sep="\t",
    header=0
)

BAM_INPUT_PATHS = MISSING_FASTQ.loc[MISSING_FASTQ["input_path"].str.endswith(".bam"), "path_hash"].values
CONSTRAINT_BAM_INPUT = "(" + "|".join(BAM_INPUT_PATHS) + ")"

FASTQ_INPUT_PATHS = MISSING_FASTQ.loc[MISSING_FASTQ["input_path"].str.endswith(".fastq"), "path_hash"].values
CONSTRAINT_FASTQ_INPUT = "(" + "|".join(FASTQ_INPUT_PATHS) + ")"


rule dump_bam_to_fastq:
    input:
        bam = lambda wildcards:
            MISSING_FASTQ.loc[
                MISSING_FASTQ["path_hash"] == wildcards.path_hash, "input_path"
            ].values[0]
    output:
        check = "check_files/{path_hash}.dump.ok"
    log:
        samtools = "log/{path_hash}.samtools.log"
    wildcard_constraints:
        path_hash = CONSTRAINT_BAM_INPUT
    conda:
        "../envs/samtools.yaml"
    threads: 6
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt,
        time_hrs = lambda wildcards, attempt: 48 * attempt
    params:
        fastq = lambda wildcards:
            MISSING_FASTQ.loc[
                MISSING_FASTQ["path_hash"] == wildcards.path_hash, "output_path"
            ].values[0]
    shell:
        "samtools bam2fq -n --threads {threads} {input.bam} 2> {log.samtools} |"
        " pigz -9 -c -p {threads} > {params.fastq} "
        " && "
        " touch {output}"


rule gzip_plain_fastq:
    input:
        fastq = lambda wildcards:
            MISSING_FASTQ.loc[
                MISSING_FASTQ["path_hash"] == wildcards.path_hash, "input_path"
            ].values[0]
    output:
        check = "check_files/{path_hash}.dump.ok"
    wildcard_constraints:
        path_hash = CONSTRAINT_FASTQ_INPUT
    conda:
        "../envs/samtools.yaml"
    threads: 6
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt,
        time_hrs = lambda wildcards, attempt: 1 * attempt
    params:
        fastq = lambda wildcards:
            MISSING_FASTQ.loc[
                MISSING_FASTQ["path_hash"] == wildcards.path_hash, "output_path"
            ].values[0]
    shell:
        "pigz -9 -c -p {threads} {input.fastq} > {params.fastq} "
        " && "
        " touch {output}"


rule run_dump_all:
    input:
        bam_to_fastq = expand(
            "check_files/{path_hash}.dump.ok",
            path_hash=BAM_INPUT_PATHS
        ),
        plain_to_fastq = expand(
            "check_files/{path_hash}.dump.ok",
            path_hash=FASTQ_INPUT_PATHS
        )
