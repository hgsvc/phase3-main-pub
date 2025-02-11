
assert SAMPLE_SEX is not None, "Need sample sex info for tig alignment"

rule hmmer_motif_search_extracted_seq:
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
        fasta = rules.fetch_tigs_from_sequence_files.output.fasta,
        motif = DIR_GLOBAL_REF.joinpath("{motif}.fasta")
    output:
        txt = DIR_PROC.joinpath(
            "extract_chrom", "hmmer_plain",
            "{sample}.hmmer.wd",
            "{sample}.{chrom}.{motif}.hmmer-out.txt"
        ),
        table = DIR_PROC.joinpath(
            "extract_chrom", "hmmer_plain",
            "{sample}.hmmer.wd",
            "{sample}.{chrom}.{motif}.hmmer-table.txt"
        ),
    benchmark:
        DIR_RSRC.joinpath(
            "extract_chrom", "hmmer_plain",
            "{sample}.{chrom}.{motif}.hmmer.rsrc"
        )
    conda:
        DIR_ENVS.joinpath("hmmer.yaml")
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
        "{input.motif} {input.fasta}"


rule run_all_hmmer_extract_chrom:
    input:
        tables = expand(
            rules.hmmer_motif_search_extracted_seq.output.table,
            sample=MALE_SAMPLES,
            motif=CHROM_Y_ALL_MOTIFS,
            chrom=["chrY"]
        )
