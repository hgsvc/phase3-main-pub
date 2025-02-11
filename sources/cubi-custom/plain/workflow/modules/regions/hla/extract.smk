
rule build_hla_cut_table:
    input:
        norm_paf = expand(
                WORKDIR_EVAL.joinpath(
                    "results/alignments/contig_to_ref/t2tv2/table",
                    "{sample}.asm-{asm_unit}.t2tv2.norm-paf.tsv.gz"
                ),
                sample=SAMPLES,
                asm_unit=MAIN_ASSEMBLY_UNITS
            )
    output:
        cut_table = DIR_RES.joinpath(
            "regions", "hla", "assembly_cut_table.{assembler}.tsv"
        )
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt
    params:
        hla_start = int(28e6),
        hla_end = int(34e6)
    run:
        import pandas as pd
        import pathlib as pl
        cut_table = []

        for paf in input.norm_paf:
            source_file = pl.Path(paf).name
            source_asm_unit = source_file.rsplit(".", 4)[0]
            sample = source_file.split(".asm-")[0]

            df = pd.read_csv(paf, sep="\t")
            df = df.loc[(df["target_name"] == "chr6") & (df["tp_align_type"] != 2), :]
            # contigs/alignments reaching less than 100 kbp into the window
            # are just ignored
            select_reach_left = (df["target_end"] - int(1e5)) > params.hla_start
            select_reach_right = (df["target_start"] - int(1e5) < params.hla_end)

            df = df.loc[select_reach_left & select_reach_right, :]
            for query, alignments in df.groupby("query_name"):
                min_q = alignments["query_start"].min()
                max_q = alignments["query_end"].max()
                min_t = alignments["target_start"].min()
                max_t = alignments["target_end"].max()
                # determine alignment orientation by majority
                align_orient = alignments.groupby("align_orient")["align_matching"].sum()
                align_orient = align_orient.idxmax()

                query_length = alignments["query_length"].values[0]
                if query_length < int(1e6):
                    # use small query sequences as is, likely garbage anyway
                    cut_table.append(
                        (sample, query, min_q, max_q, align_orient, min_t, max_t,
                         0, query_length, query_length, source_asm_unit, source_file))
                    continue

                # following: slack of 1 Mbp because rare
                # cases exist where the alignment starts a few
                # 100 kbp into the MHC window; assert that is does
                # not start/end in the middle of the MHC locus;
                # if that happens, need to investigate
                assert params.hla_start + int(1e6) > min_t, f"{sample} / {query} / {min_t}"
                offset_start = params.hla_start - min_t
                assert params.hla_end - int(1e6) < max_t, f"{sample} / {query} / {max_t}"
                offset_end = params.hla_end - max_t

                if align_orient > 0:
                    # Alignment in forward orientation;
                    # NB: offset end will typically be negative
                    cut_begin = min_q + offset_start
                    cut_end = max_q + offset_end
                elif align_orient < 0:
                    # Alignment in reverse orientation,
                    # hence end in query is start in target
                    # and vice versa, so flip sign of offset
                    # values
                    cut_begin = min_q + offset_end * -1
                    cut_end = max_q + offset_start * -1
                else:
                    raise ValueError("Cannot process alignment orientation 0")
                cut_length = cut_end - cut_begin
                assert cut_length > int(5.5e6)
                assert cut_length < int(6.5e6)
                cut_table.append(
                    (sample, query, min_q, max_q, align_orient, min_t, max_t,
                     cut_begin, cut_end, cut_length,
                     source_asm_unit, source_file)
                )

        cut_table = pd.DataFrame.from_records(
            cut_table, columns=[
                "sample", "query", "query_start", "query_end",
                "query_orientation", "target_start", "target_end",
                "cut_query_begin", "cut_query_end",
                "cut_length", "source_assembly", "source_file"
            ]
        )
        cut_table.to_csv(output.cut_table, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule extract_hla_sequences:
    input:
        fasta = get_asm_unit_fasta_files,
        cut_table = rules.build_hla_cut_table.output.cut_table
    output:
        table_comp = DIR_RES.joinpath(
            "regions", "hla", "assembly_seq_composition.{assembler}.tsv"
        ),
        all_seqs = DIR_RES.joinpath(
            "regions", "hla", "assembly_all_hla.{assembler}.fasta.gz"
        )
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 4096 * attempt
    params:
        script=DIR_SCRIPTS.joinpath("extract_region.py"),
        out_dir=lambda wildcards, output: pathlib.Path(output.table_comp).parent.joinpath("single_fasta")
    shell:
        "{params.script} --fasta {input.fasta} --cut-table {input.cut_table} "
            "--add-suffix HLA --dump-separately {params.out_dir} --clean-out-dir "
            "--dump-merged stdout --dump-stats {output.table_comp} "
            " | "
        "bgzip > {output.all_seqs}"
            " && "
        "samtools faidx {output.all_seqs}"


rule extract_hla_reference_sequence:
    input:
        fasta = WORKDIR_EVAL.joinpath(
            "global_ref", "chm13v2.0.fasta"
        ),
        faidx = WORKDIR_EVAL.joinpath(
            "global_ref", "chm13v2.0.fasta.fai"
        )
    output:
        fasta = DIR_LOCAL_REF.joinpath(
            "chm13v2.0.HLA.fasta.gz"
        )
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    params:
        hla_start = int(28e6),
        hla_end = int(34e6)
    shell:
        "samtools faidx {input.fasta} chr6:{params.hla_start}-{params.hla_end} "
        " | gzip > {output.fasta}"


rule blast_check_extracted_sequences:
    """NB: blast cannot read gzipped files ...
    """
    input:
        all_seqs = rules.extract_hla_sequences.output.all_seqs,
        ref_seq = rules.extract_hla_reference_sequence.output.fasta,
        hla = DIR_GLOBAL_REF.joinpath("hla_coding_transcripts.nuc.fasta")
    output:
        blast_all = DIR_PROC.joinpath("regions", "hla", "blast_out", "assembly_all_hla.blast.{assembler}.txt"),
        blast_ref = DIR_PROC.joinpath("regions", "hla", "blast_out", "chm13v2.0_hla.blast.{assembler}.txt"),
    conda:
        DIR_ENVS.joinpath("seqtools.yaml")
    params:
        tmp_all=lambda wildcards, input: pathlib.Path(input.all_seqs).with_suffix(".tmp.fa"),
        tmp_ref=lambda wildcards, input: pathlib.Path(input.ref_seq).with_suffix(".tmp.fa"),
        word_size=101,
        evalue=0.001
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt,
    shell:
        "pigz -d -c {input.all_seqs} > {params.tmp_all}"
            " && "
        "pigz -d -c {input.ref_seq} > {params.tmp_ref}"
            " && "
        "blastn -query {input.hla} -subject {params.tmp_ref} "
            "-word_size {params.word_size} -outfmt=6 -evalue {params.evalue} "
        "-out {output.blast_ref} && rm {params.tmp_ref} "
        " && "
        "blastn -query {input.hla} -subject {params.tmp_all} "
            "-word_size {params.word_size} -outfmt=6 -evalue {params.evalue} "
        "-out {output.blast_all} && rm {params.tmp_all}"


localrules: summarize_blast_hits
rule summarize_blast_hits:
    input:
        blast_all = rules.blast_check_extracted_sequences.output.blast_all,
        blast_ref = rules.blast_check_extracted_sequences.output.blast_ref
    output:
        summary = DIR_RES.joinpath(
            "regions", "hla", "blast_hits.{assembler}.summary.tsv.gz"
        ),
        table = DIR_RES.joinpath(
            "regions", "hla", "blast_hits.{assembler}.tsv.gz"
        )
    run:
        import pandas as pd
        blast_columns = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore".split()
        df = pd.concat(
            [
                pd.read_csv(input.blast_ref, sep="\t", header=None, names=blast_columns),
                pd.read_csv(input.blast_all, sep="\t", header=None, names=blast_columns)
            ],
            axis=0, ignore_index=False
        )
        df.sort_values(["sseqid", "bitscore"], ascending=[True, False], inplace=True)
        df.to_csv(output.table, sep="\t", header=True, index=False)

        summary = df.groupby(["sseqid", "qseqid"]).agg(
            median_bitscore=("bitscore", "median"),
            median_evalue=("evalue", "median"),
            total_hits=("length", "count")
        )
        summary.sort_index(inplace=True)
        summary.to_csv(output.summary, sep="\t", header=True, index=True)
    # END OF RUN BLOCK


localrules: blast_summarize_hits_per_sequence
rule blast_summarize_hits_per_sequence:
    input:
        summary = rules.summarize_blast_hits.output.summary,
        all_seqs = rules.extract_hla_sequences.output.all_seqs
    output:
        hit_count = DIR_RES.joinpath(
            "regions", "hla", "blast_hits_per_seq.{assembler}.tsv"
        )
    run:
        import pathlib as pl
        import pandas as pd
        summary = pd.read_csv(input.summary, sep="\t")
        summary["sample"] = summary["sseqid"].apply(lambda x: "chm13v2.0" if "chr6" in x else x.split("_")[0])

        def assign_hap(seqid):
            if "haplotype1" in seqid:
                return "asm-hap1"
            if "haplotype2" in seqid:
                return "asm-hap2"
            if "unassigned" in seqid:
                return "asm-unassigned"
            return "chm13v2.0"
        summary["assembly"] = summary["sseqid"].apply(assign_hap)

        single_seq_folder = pl.Path(input.all_seqs).parent.joinpath("single_fasta")
        single_seq_fastas = sorted([filepath.name for filepath in single_seq_folder.glob("*.fasta.gz")])
        fasta_df = []
        for fasta in single_seq_fastas:
            sample, asm_unit, _, _, _ = fasta.rsplit(".", 4)
            fasta_df.append((sample, asm_unit, fasta))
        fasta_df.append(("chm13v2.0", "chm13v2.0", "chm13v2.0.HLA.fasta.gz"))
        fasta_df = pd.DataFrame.from_records(fasta_df, columns=["sample", "assembly", "fasta"])
        summary = summary.merge(fasta_df, on=["sample", "assembly"], how="outer")
        summary["total_hits"].fillna(0, inplace=True)
        summary.rename({"total_hits": "total_blast_hla_hits"}, axis=1, inplace=True)
        hits_per_seq = summary.groupby(["sample", "assembly", "fasta"])["total_blast_hla_hits"].sum()
        hits_per_seq.sort_index(inplace=True)
        hits_per_seq.to_csv(output.hit_count, sep="\t", header=True, index=True)
    # END OF RUN BLOCK


# Irrelevant rule - just kept for reference
# rule produce_hla_msa:
#     input:
#         ref_fasta = rules.extract_hla_reference_sequence.output.fasta,
#         asm_fasta = rules.extract_hla_sequences.output.all_seqs,
#     output:
#         msa = DIR_RES.joinpath(
#             "regions", "hla", "assembly_all_hla.msa.fasta"
#         ),
#         distmat = DIR_RES.joinpath(
#             "regions", "hla", "assembly_all_hla.msa-distmat.txt"
#         )
#     log:
#         DIR_LOG.joinpath("regions", "hla", "extract", "clustalo.log")
#     conda:
#         DIR_ENVS.joinpath("msa.yaml")
#     threads: 12
#     resources:
#         mem_mb = lambda wildcards, attempt: 24576 * attempt,
#         time_hrs = lambda wildcards, attempt: 47 * attempt
#     shell:
#         "zcat {input} | clustalo --in - --seqtype=DNA --infmt=fasta "
#         "--percent-id --outfmt=fasta --threads {threads} --verbose "
#         "--full --distmat-out={output.distmat} "
#         "--out {output.msa} &> {log}"


rule run_blast_checks:
    input:
        summary = rules.blast_summarize_hits_per_sequence.output.hit_count


rule copy_files_to_share:
    input:
        cut_table = rules.build_hla_cut_table.output.cut_table,
        all_seqs = rules.extract_hla_sequences.output.all_seqs,
        seq_comp = rules.extract_hla_sequences.output.table_comp,
        ref_seq = rules.extract_hla_reference_sequence.output.fasta,
        blast_summary = rules.summarize_blast_hits.output.summary,
        blast_table = rules.summarize_blast_hits.output.table,
        blast_count = rules.blast_summarize_hits_per_sequence.output.hit_count
    output:
        file_listing = DIR_PROC.joinpath(
            "regions", "hla", "files_shared.{assembler}.lst"
        )
    params:
        single_fasta=lambda wildcards, input: pathlib.Path(input.seq_comp).parent.joinpath("single_fasta")
    run:
        import pathlib as pl
        import shutil as sh
        HLA_GLOBUS_SHARE = BASE_SHARE_LOCATION.joinpath("sig_hla")
        HLA_GLOBUS_SHARE.mkdir(exist_ok=True, parents=True)
        index_extensions = [".fai", ".tbi", ".gzi"]

        copied_files = []
        for share_file in input:
            source_path = pl.Path(share_file)
            assert source_path.is_file()
            target_path = HLA_GLOBUS_SHARE.joinpath(source_path.name)
            sh.copy(source_path, target_path)
            copied_files.append(
                (str(source_path), str(target_path))
            )
            for index_ext in index_extensions:
                source_ext = source_path.suffix  # starts with "."
                source_index = source_path.with_suffix(f"{source_ext}{index_ext}")
                if source_index.is_file():
                    target_index = HLA_GLOBUS_SHARE.joinpath(source_index.name)
                    sh.copy(source_index, target_index)
                    copied_files.append(
                        (str(source_index), str(target_index))
                    )

        target_single_seq = HLA_GLOBUS_SHARE.joinpath("single_fasta")
        target_single_seq.mkdir(exist_ok=True, parents=True)
        source_single_seq = pl.Path(params.single_fasta)
        for source_file in source_single_seq.iterdir():
            if not source_file.is_file():
                continue
            target_file = target_single_seq.joinpath(source_file.name)
            sh.copy(source_file, target_file)

        with open(output.file_listing, "w") as listing:
            for src, trg in copied_files:
                listing.write(f"{src}\t{trg}\n")
    # END OF RUN BLOCK


rule run_all_extract_hla:
    input:
        seqs = expand(
            rules.extract_hla_sequences.output.all_seqs,
            assembler=ASSEMBLER
        ),
        blast = expand(
            rules.blast_summarize_hits_per_sequence.output.hit_count,
            assembler=ASSEMBLER
        ),
        file_listing = expand(
            rules.copy_files_to_share.output.file_listing,
            assembler=ASSEMBLER
        )
