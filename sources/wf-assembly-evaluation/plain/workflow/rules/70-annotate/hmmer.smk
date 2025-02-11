
rule hmmer_motif_search:
    """NB: the reported hits can EITHER
    be thresholded on the E-value [-E] OR
    on the score [-T], but not both in the same run.
    The implementation here only thresholds on
    the E-value (if specified for the motif) and
    then later labels hits above the score threshold
    (if specified for the motif) as high-quality

    IMPORTANT:
    Only v3.4+ of HMMER has a bug fix for an invalid
    alphabet detection:
    https://github.com/EddyRivasLab/hmmer/pull/252
    """
    input:
        asm_unit = rules.create_plain_assembly_file.output.tmp_fa,
        motif = DIR_GLOBAL_REF.joinpath("{motif}.fasta")
    output:
        txt = DIR_PROC.joinpath(
            "70-annotate", "hmmer",
            "{sample}.{asm_unit}.hmmer.wd",
            "{sample}.{asm_unit}.{motif}.hmmer-out.txt"
        ),
        table = DIR_PROC.joinpath(
            "70-annotate", "hmmer",
            "{sample}.{asm_unit}.hmmer.wd",
            "{sample}.{asm_unit}.{motif}.hmmer-table.txt"
        ),
    log:
        DIR_LOG.joinpath(
            "70-annotate", "hmmer",
            "{sample}.{asm_unit}.{motif}.hmmer.log"
        )
    benchmark:
        DIR_RSRC.joinpath(
            "70-annotate", "hmmer",
            "{sample}.{asm_unit}.{motif}.hmmer.rsrc"
        )
    conda:
        DIR_ENVS.joinpath("biotools", "hmmer.yaml")
    threads: lambda wildcards: get_num_threads_hmmer(wildcards.motif)
    resources:
        mem_mb = lambda wildcards, attempt: attempt * get_mem_mb_hmmer(wildcards.motif),
        time_hrs = lambda wildcards, attempt: attempt * hmmer_scaling("time", wildcards.motif)
    params:
        evalue_t = lambda wildcards: hmmer_threshold_value("evalue_t", wildcards.motif),
    shell:
        "nhmmer --cpu {threads} --dna "
        "-o {output.txt} --tblout {output.table} "
        "-E {params.evalue_t} "
        "{input.motif} {input.asm_unit} &> {log}"


rule normalize_hmmer_table_output:
    input:
        table = rules.hmmer_motif_search.output.table,
        motif = DIR_GLOBAL_REF.joinpath("{motif}.fasta")
    output:
        norm = DIR_RES.joinpath(
            "annotations", "hmmer",
            "{sample}.{asm_unit}.{motif}.hmmer.matches.tsv.gz"

        ),
        agg = DIR_RES.joinpath(
            "annotations", "hmmer",
            "{sample}.{asm_unit}.{motif}.hmmer.agg-hits.tsv.gz"

        ),
        bed = DIR_PROC.joinpath(
            "70-annotate", "hmmer",
            "{sample}.{asm_unit}.{motif}.hmmer.hits-all.bed"

        ),
        hiq = DIR_PROC.joinpath(
            "70-annotate", "hmmer",
            "{sample}.{asm_unit}.{motif}.hmmer.hits-hiq.bed"
        )
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    params:
        script=find_script("norm_hmmer"),
        score_t = lambda wildcards: hmmer_threshold_value("score_t", wildcards.motif),
        acc_res=lambda wildcards, output: register_result(output.norm, output.agg)
    shell:
        "{params.script} --hmmer-table {input.table} --motif-file {input.motif} "
            "--score-threshold {params.score_t} --force-empty-output --add-empty-flag "
            "--norm-table {output.norm} --aggregated {output.agg} "
            "--bed-all {output.bed} --bed-hiq {output.hiq}"


rule compress_index_hmmer_bed_hits:
    """The tabix indexing of the bed files
    would fail on empty input. This rule
    checks if the .EMPTY flag file does not exist
    and if so, proceeds with bgzip and tabix

    bash: ! -f => EMPTY flag file does not exist
    """
    input:
        bed = DIR_PROC.joinpath(
            "70-annotate", "hmmer",
            "{sample}.{asm_unit}.{motif}.hmmer.{bed_type}.bed"
        )
    output:
        bgz = DIR_RES.joinpath(
            "annotations", "hmmer",
            "{sample}.{asm_unit}.{motif}.hmmer.{bed_type}.bed.gz"
        ),
        tbi = DIR_RES.joinpath(
            "annotations", "hmmer",
            "{sample}.{asm_unit}.{motif}.hmmer.{bed_type}.bed.gz.tbi"
        ),
    wildcard_constraints:
        bed_type="(hits\-all|hits\-hiq)"
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    params:
        acc_res=lambda wildcards, output: register_result(output)
    shell:
        "rm -f {output.bgz}.EMPTY ; rm -f {output.tbi}.EMPTY ; "
        "if [ ! -f {input.bed}.EMPTY ] ; then "
        "{{ "
        "cat {input.bed} | bgzip > {output.bgz}"
            " && "
        "tabix -p bed {output.bgz} ; "
        "}} else {{ "
        "touch {output.bgz} && touch {output.tbi} "
        "&& touch {output.bgz}.EMPTY && touch {output.tbi}.EMPTY ; "
        "}} fi;"


rule run_all_hmmer_jobs:
    input:
        tables = expand(
            rules.normalize_hmmer_table_output.output.norm,
            sample=SAMPLES,
            asm_unit=ASSEMBLY_UNITS_NO_CONTAM,
            motif=HMMER_MOTIF_NAMES
        ),
        stats = expand(
            rules.normalize_hmmer_table_output.output.agg,
            sample=SAMPLES,
            asm_unit=ASSEMBLY_UNITS_NO_CONTAM,
            motif=HMMER_MOTIF_NAMES
        ),
        bed_out = expand(
            rules.compress_index_hmmer_bed_hits.output.bgz,
            sample=SAMPLES,
            asm_unit=ASSEMBLY_UNITS_NO_CONTAM,
            motif=HMMER_MOTIF_NAMES,
            bed_type=["hits-all", "hits-hiq"]
        ),
