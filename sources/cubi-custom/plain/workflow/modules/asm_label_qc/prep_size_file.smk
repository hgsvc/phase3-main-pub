
if ASSEMBLER == "verkko":

    localrules: prep_verkko_fasta_index
    rule prep_verkko_fasta_index:
        input:
            fai = expand(
                WORKDIR_EVAL.joinpath(
                    "results/assemblies/{sample}",
                    "{sample}.asm-{asm_unit}.fasta.gz.fai"
                ),
                asm_unit=MAIN_ASSEMBLY_UNITS,
                allow_missing=True
            )
        output:
            sizes = DIR_PROC.joinpath(
                "asm_label_qc", "assembly_size",
                "{sample}.sizes.txt"
            )
        run:
            import io
            out = io.StringIO()
            for fai_file in input.fai:
                with open(fai_file) as listing:
                    for line in listing:
                        seq, size = line.split()[:2]
                        _ = int(size)
                        out.write(f"{seq}\t{size}\n")
            with open(output.sizes, "w") as dump:
                _ = dump.write(out.getvalue())
        # END OF RUN BLOCK


if ASSEMBLER == "hifiasm":

    localrules: prep_hifiasm_fasta_index
    rule prep_hifiasm_fasta_index:
        input:
            fai = WORKDIR_EVAL.joinpath(
                "proc/05-preprocess/merge_tag_asm",
                "{sample}.asm-mrg-tag.fasta.fai"
            )
        output:
            sizes = DIR_PROC.joinpath(
                "asm_label_qc", "assembly_size",
                "{sample}.sizes.txt"
            )
        run:
            import io
            out = io.StringIO()
            with open(input.fai) as listing:
                for line in listing:
                    seq, size = line.split()[:2]
                    # sequence is suffixed with tag as
                    # part of preprocessing
                    seq = seq.rsplit(".", 1)[0]
                    _ = int(size)
                    out.write(f"{seq}\t{size}\n")
            with open(output.sizes, "w") as dump:
                _ = dump.write(out.getvalue())
        # END OF RUN BLOCK


rule run_all_prep_size_files:
    input:
        listings = expand(
            DIR_PROC.joinpath(
                "asm_label_qc", "assembly_size",
                "{sample}.sizes.txt"
            ),
            sample=SAMPLES
        )
