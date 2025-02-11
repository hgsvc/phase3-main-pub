
rule rename_extracted_tigs:
    input:
        fasta = rules.fetch_tigs_from_sequence_files.output.fasta,
        sqcls_align = DIR_RES.joinpath(
            "extract_chrom", "seqclass_to_chrom",
            "minimap", "{sample}.t2t.{chrom}.cls.norm-paf.tsv.gz"
        ),
        ref_align = DIR_RES.joinpath(
            "extract_chrom", "chrom_to_ref",
            "minimap", "{sample}.t2t.{chrom}.seq.norm-paf.tsv.gz"
        ),
        sqcls_order = DIR_GLOBAL_REF.joinpath("T2T.chrY-seq-classes.tsv")
    output:
        fasta = DIR_RES.joinpath(
            "extract_chrom", "renamed",
            "fasta_seqs", "{sample}.{chrom}.renamed.fasta"
        ),
        sqcls_align = DIR_RES.joinpath(
            "extract_chrom", "renamed", "seqclass_to_chrom",
            "minimap", "{sample}.t2t.{chrom}.cls.norm-paf.tsv.gz"
        ),
        ref_align = DIR_RES.joinpath(
            "extract_chrom", "renamed", "chrom_to_ref",
            "minimap", "{sample}.t2t.{chrom}.seq.norm-paf.tsv.gz"
        ),
        table = DIR_RES.joinpath(
            "extract_chrom", "renamed",
            "{sample}.t2t.{chrom}.rename-info.tsv"
        )
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    params:
        script = DIR_SCRIPTS.joinpath("extract_chrom", "rename_tigs.py")
    shell:
        "{params.script} --min-align-total 25000 --input-class-order {input.sqcls_order} "
        "--fasta-in {input.fasta} --fasta-out {output.fasta} "
        "--class-align {input.sqcls_align} --out-class-align {output.sqcls_align} "
        "--seq-align {input.ref_align} --out-seq-align {output.ref_align} "
        "--out-table {output.table}"


localrules: rename_dump_sed_file
rule rename_dump_sed_file:
    input:
        table = rules.rename_extracted_tigs.output.table
    output:
        txt = DIR_RES.joinpath(
            "extract_chrom", "renamed", "aux",
            "{sample}.t2t.{chrom}.ren-old-new.sed"
        )
    run:
        import pandas as pd
        df = pd.read_csv(input.table, sep="\t", header=0, comment="#")
        with open(output.txt, "w") as dump:
            for row in df.itertuples():
                dump.write(f"s/\\b{row.old_name}\\b/{row.new_name}/g\n")
    # END OF RUN BLOCK


rule rename_hmmer_motif_hits:
    input:
        txt = rules.hmmer_motif_search_extracted_seq.output.txt,
        table = rules.hmmer_motif_search_extracted_seq.output.table,
        sed = rules.rename_dump_sed_file.output.txt
    output:
        txt = DIR_RES.joinpath(
            "extract_chrom", "renamed", "motif_hits",
            "{sample}", "{sample}.{chrom}.{motif}.hmmer-out.txt.gz"
        ),
        table = DIR_RES.joinpath(
            "extract_chrom", "renamed", "motif_hits",
            "{sample}", "{sample}.{chrom}.{motif}.hmmer-table.txt.gz"
        ),
    shell:
        "sed -f {input.sed} {input.txt} | gzip > {output.txt}"
            " && "
        "sed -f {input.sed} {input.table} | gzip > {output.table}"


rule run_all_rename_extracted_chrom:
    input:
        tables = expand(
            rules.rename_extracted_tigs.output.table,
            sample=MALE_SAMPLES,
            chrom=["chrY"]
        ),
        sed_txt = expand(
            rules.rename_dump_sed_file.output.txt,
            sample=MALE_SAMPLES,
            chrom=["chrY"]
        ),
        hmmer = expand(
            rules.rename_hmmer_motif_hits.output.table,
            sample=MALE_SAMPLES,
            chrom=["chrY"],
            motif=CHROM_Y_ALL_MOTIFS
        )
