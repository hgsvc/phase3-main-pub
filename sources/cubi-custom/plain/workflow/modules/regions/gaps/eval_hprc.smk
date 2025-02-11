######################################
# NB - This module is outdated
# and has been replaced by
# 'annotate_gaps.smk', which includes
# the Flagger and NucFreq error
# annotation
# JUST HERE FOR REFERENCE
######################################


import pathlib

# HPRC gap file must exist
# NB: is in chm13/T2T coordinates
HPRC_COMMON_GAPS = DIR_ANNOTATIONS.joinpath(
    "roi", "porubsky2023_hprc_common-gaps.tsv"
).resolve(strict=True)

# some sample karyotype assignment must exist
ASM_UNIT_SEX_ALL = WORKDIR_EVAL.joinpath(
    "results/reports/ref_chrom_assign",
    "karyo-est.hgsvc3-and-special.tsv"
)

ASM_UNIT_SEX_HGSVC3 = WORKDIR_EVAL.joinpath(
    "results/reports/ref_chrom_assign",
    "karyo-est.hgsvc3.tsv"
)

if ASM_UNIT_SEX_HGSVC3.is_file():
    ASM_UNIT_SEX_FILE = ASM_UNIT_SEX_HGSVC3
elif ASM_UNIT_SEX_ALL.is_file():
    ASM_UNIT_SEX_FILE = ASM_UNIT_SEX_ALL
else:
    raise RuntimeError("Need AU sex annotation")


ASM_UNIT_SEX = load_assembly_unit_karyotypes(ASM_UNIT_SEX_FILE)


rule annotate_gaps_contig_cov:
    input:
        gaps = HPRC_COMMON_GAPS,
        ctg_cov = WORKDIR_EVAL.joinpath(
            "results/coverage/contig_ref/t2tv2",
            "{sample}.t2tv2.win-ctg-cov.tsv.gz"
        )
    output:
        agg = DIR_RES.joinpath(
            "regions", "gaps", "hprc",
            "sample_ctg_cov", "{sample}.ctg-cov.hprc-common.tsv"
        )
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    resources:
        time_hrs=0
    params:
        script=DIR_SCRIPTS.joinpath("eval_gaps", "annotate_ctg_cov.py"),
        sex_label=lambda wildcards: ASM_UNIT_SEX[wildcards.sample]
    shell:
        "{params.script} --ctg-cov {input.ctg_cov} --gaps {input.gaps} "
        "--out-agg {output.agg} --sample {wildcards.sample} --assembly-sex {params.sex_label}"


rule label_gaps_contig_cov:
    input:
        table = rules.annotate_gaps_contig_cov.output.agg
    output:
        table = DIR_RES.joinpath(
            "regions", "gaps", "hprc",
            "labeled_ctg_cov", "{sample}.ctg-cov.hprc-common.labeled.tsv"
        )
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    resources:
        time_hrs=0
    params:
        script=DIR_SCRIPTS.joinpath("eval_gaps", "label_gap_cov.py"),
        aln_label="bplvl"
    shell:
        "{params.script} --input-table {input.table} --output-table {output.table} "
        "--aln-label {params.aln_label}"


rule run_all_hprc_gaps:
    input:
        ctg_cov_label = expand(
            rules.label_gaps_contig_cov.output.table,
            sample=SAMPLES
        )
