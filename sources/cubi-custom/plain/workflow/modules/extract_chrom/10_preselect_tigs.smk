
rule preselect_chrom_contigs:
    input:
        aln = expand(
            WORKDIR_EVAL.joinpath(
                "results/alignments/contig_to_ref/t2tv2/table",
                "{{sample}}.asm-{asm_unit}.t2tv2.norm-paf.tsv.gz"
            ),
            asm_unit=MAIN_ASSEMBLY_UNITS
        ),
        motif_hits = expand(
            WORKDIR_EVAL.joinpath(
                "results", "annotations", "hmmer",
                "{{sample}}.asm-{asm_unit}.{motif}.hmmer.agg-hits.tsv.gz"
            ),
            asm_unit=MAIN_ASSEMBLY_UNITS,
            motif=CHROM_Y_SELECT_MOTIFS
        )
    output:
        contig_list = DIR_RES.joinpath(
            "extract_chrom", "select_tigs",
            "{sample}.{chrom}.selected_tigs.tsv"
        )
    wildcard_constraints:
        chrom = "(chrX|chrY)"
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    params:
        script = DIR_SCRIPTS.joinpath("extract_chrom", "preselect_contigs.py"),
        filter_aln = 2,
        min_aln_t = 50
    shell:
        "{params.script} --contig-ref-align {input.aln} --motif-hits {input.motif_hits} "
            "--select-chrom {wildcards.chrom} --drop-alignments {params.filter_aln} "
            "--min-aln-threshold {params.min_aln_t} --output {output.contig_list}"


rule fetch_tigs_from_sequence_files:
    input:
        listing = rules.preselect_chrom_contigs.output.contig_list,
        fastas = get_asm_unit_fasta_files
    output:
        fasta = DIR_PROC.joinpath(
            "extract_chrom", "fasta_seqs",
            "{sample}.{chrom}.oriented.fasta"
        )
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    wildcard_constraints:
        chrom="(chrY|chrX)"
    params:
        script=DIR_SCRIPTS.joinpath("extract_chrom", "fetch_seq.py")
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    shell:
        "{params.script} --fasta-files {input.fastas} --selected-contigs {input.listing} "
            "--output {output.fasta}"


rule run_all_fetch_chrom_contigs:
    input:
        fasta = expand(
            rules.fetch_tigs_from_sequence_files.output.fasta,
            sample=SAMPLES,
            chrom=["chrY", "chrX"]
        )
