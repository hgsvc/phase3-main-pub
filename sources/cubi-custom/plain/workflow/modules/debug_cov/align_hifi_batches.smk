"""
This module is only executed on demand
to check for uneven chromosome coverage
across data batches, which could be
indicative of a cell line artifact
"""

import pathlib as pl

DATA_ROOT = pl.Path("/gpfs/project/projects/medbioinf/data/00_RESTRUCTURE/project-centric/hgsvc/pacbio_hifi")

# so far, all female
ALIGN_SAMPLES = [
    "HG00732",
    "HG00514",
    "NA18939"
]

BATCHES_PER_SAMPLE = {
    "HG00732": {
        "batch1_2019": [
            "20190925_PUR_PacBio_HiFi/HG00732_20190925_EEE_m54329U_190604_224858.Q20.fastq.gz",
            "20190925_PUR_PacBio_HiFi/HG00732_20190925_EEE_m54329U_190610_071123.Q20.fastq.gz",
            "20190925_PUR_PacBio_HiFi/HG00732_20190925_EEE_m54329U_190611_132246.Q20.fastq.gz",
            "20190925_PUR_PacBio_HiFi/HG00732_20190925_EEE_m54329U_190612_193411.Q20.fastq.gz",
            "20190925_PUR_PacBio_HiFi/HG00732_20190925_EEE_m54329U_190705_000551.Q20.fastq.gz",
            "20200722_PUR_PacBio_HG00732_HiFi/HG00732_20200722_EEE_m54329U_200528_200534.ccs.fastq.gz",
            "20200722_PUR_PacBio_HG00732_HiFi/HG00732_20200722_EEE_m64076_200601_234627.ccs.fastq.gz",
            "20200722_PUR_PacBio_HG00732_HiFi/HG00732_20200722_EEE_m64076_200603_055852.ccs.fastq.gz"
        ],
        "batch2_2020": [
            "20200722_PUR_PacBio_HG00732_HiFi/HG00732_20200722_EEE_m54329U_200528_200534.ccs.fastq.gz",
            "20200722_PUR_PacBio_HG00732_HiFi/HG00732_20200722_EEE_m64076_200601_234627.ccs.fastq.gz",
            "20200722_PUR_PacBio_HG00732_HiFi/HG00732_20200722_EEE_m64076_200603_055852.ccs.fastq.gz"
        ],
        "batch3_2023": [
            "20231126_UW_HiFi/HG00732/m84046_231113_224421_s1.hifi_reads.bc2033.fastq.gz",
            "20231126_UW_HiFi/HG00732/m84046_231117_203108_s3.hifi_reads.bc2033.fastq.gz"
        ]
    },
    "HG00514": {
        "batch1_2020": [
            "20200731_CHS_PacBio_HG00514_HiFi_reseq/m54329U_200715_194535.ccs.fastq.gz",
            "20200731_CHS_PacBio_HG00514_HiFi_reseq/m54329U_200717_235548.ccs.fastq.gz",
            "20200731_CHS_PacBio_HG00514_HiFi_reseq/m54329U_200719_061020.ccs.fastq.gz"
        ],
        "batch2_2023": [
            "20231126_UW_HiFi/HG00514/m84046_231007_211119_s1.hifi_reads.bc2068.fastq.gz"
        ]
    },
    "NA18939": {
        "batch1_2022": [
            "20230512_HGSVC_EEE_HIFI/NA18939/m64076_220402_104632.hifi_reads.fastq.gz",
            "20230512_HGSVC_EEE_HIFI/NA18939/m64076_220403_214315.hifi_reads.fastq.gz",
            "20230512_HGSVC_EEE_HIFI/NA18939/m64076_220405_084052.hifi_reads.fastq.gz"
        ],
        "batch2_2023": [
            "20230512_HGSVC_EEE_HIFI/NA18939/m64076_230717_214340-bc2082.hifi_reads.fastq.gz"
        ]
    }
}


rule align_batches:
    input:
        reads = lambda wildcards: [DATA_ROOT.joinpath(fp) for fp in BATCHES_PER_SAMPLE[wildcards.sample][wildcards.batch]],
        ref = DIR_GLOBAL_REF.joinpath("chm13v2.0_maskedY_rCRS.fa.gz"),
        ref_idx = DIR_GLOBAL_REF.joinpath("chm13v2.0_maskedY_rCRS.fa.gz.fai"),
        bed = DIR_GLOBAL_REF.joinpath("chm13v2.0_dip-chroms.bed")
    output:
        bam = DIR_PROC.joinpath(
            "debug_cov/alns/{sample}/{sample}_{batch}.sort.bam"
        ),
        bai = DIR_PROC.joinpath(
            "debug_cov/alns/{sample}/{sample}_{batch}.sort.bam.bai"
        )
    conda:
        DIR_ENVS.joinpath("minimap.yaml")
    threads: 24
    resources:
        time_hrs=lambda wildcards, attempt: 23 * attempt,
        mem_mb=lambda wildcards, attempt: 100352 * attempt
    shell:
        "minimap2 -x map-hifi --secondary=no -a --eqx -t {threads} {input.ref} {input.reads}"
            " | "
        "samtools view -F 1796 -L {input.bed} --with-header"
            " | "
        "samtools sort --threads 6 -l 9 -m 1G "
        "-T tmp/{wildcards.sample}_{wildcards.batch} "
        "-o {output.bam}"
            " && "
        "samtools index {output.bam}"


rule compute_coverage:
    input:
        bam = rules.align_batches.output.bam,
        bai = rules.align_batches.output.bai,
    output:
        txt = DIR_PROC.joinpath(
            "debug_cov/{sample}/{sample}_{batch}.cov.txt"
        )
    conda:
        DIR_ENVS.joinpath("minimap.yaml")
    resources:
        time_hrs=lambda wildcards, attempt: attempt,
        mem_mb=lambda wildcards, attempt: 4096 * attempt
    shell:
        "samtools coverage --depth 99 -o {output.txt} {input.bam}"


localrules: check_all_batches_complete
rule check_all_batches_complete:
    input:
        txt = lambda wildcards: expand(
            rules.compute_coverage.output.txt,
            batch=sorted(BATCHES_PER_SAMPLE[wildcards.sample].keys()),
            allow_missing=True
        )
    output:
        ok = DIR_PROC.joinpath(
            "debug_cov", "{sample}.all-batches-complete.ok"
        )
    shell:
        "touch {output.ok}"


rule run_all_debug_cov:
    input:
        cov = expand(
            rules.check_all_batches_complete.output.ok,
            sample=ALIGN_SAMPLES
        )
