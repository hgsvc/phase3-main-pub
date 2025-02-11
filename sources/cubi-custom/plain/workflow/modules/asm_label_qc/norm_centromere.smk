
localrules: split_centromere_annotation
rule split_centromere_annotation:
    input:
        listing = CENTROMERE_ANNOTATION
    output:
        beds = expand(
            DIR_RES.joinpath(
                "asm_label_qc", "norm_tables", "centromeres",
                "{sample}.active_asat_HOR_arrays_v3.bed"
            ),
            sample=SAMPLES
        )
    run:
        import pandas as pd
        header = ["sample_plain", "seq_loc", "chrom", "strand"]
        df = pd.read_csv(input.listing, sep="\t", header=None, names=header)
        def get_seq(seq_loc):
            seq, loc = seq_loc.split(":")
            start, end = loc.split("-")
            start = int(start)
            end = int(end)
            assert start < end
            return seq, start, end
        seqs = df["seq_loc"].apply(get_seq)
        seqs = pd.DataFrame.from_records(seqs, index=df.index, columns=["seq", "start", "end"])
        df = df.join(seqs)
        df.drop("seq_loc", axis=1, inplace=True)
        def assign_assembler(seq):
            if any(x in seq for x in ["h1tg", "h2tg"]):
                return ".hsm-ps-sseq"
            elif any(x in seq for x in ["haplotype1", "haplotype2", "unassigned"]):
                return ".vrk-ps-sseq"
            else:
                raise ValueError(seq)
        df["assembler"] = df["seq"].apply(assign_assembler)
        df["sample"] = df["sample_plain"] + df["assembler"]
        df = df[["seq", "start", "end", "chrom", "strand", "sample"]]
        df.insert(4, "score", 1)
        df.sort_values(["sample", "seq", "start", "end"], inplace=True)
        df.rename({"seq": "#seq"}, axis=1, inplace=True)

        ann_samples = set()
        for sample, cens in df.groupby("sample"):
            out_file = DIR_RES.joinpath(
                "asm_label_qc", "norm_tables", "centromeres",
                f"{sample}.active_asat_HOR_arrays_v3.bed"
            )
            cens.to_csv(out_file, sep="\t", header=True, index=False)
            ann_samples.add(sample)

        remain = set(SAMPLES) - ann_samples
        if remain:
            for sample in sorted(remain):
                gsize = DIR_PROC.joinpath(
                    "asm_label_qc", "assembly_size",
                    f"{sample}.sizes.txt"
                ).resolve(strict=True)
                mock = pd.read_csv(gsize, sep="\t", header=None, names=["seq", "end"])
                mock["start"] = 0
                mock["score"] = 0
                mock["chrom"] = "chrUN"
                mock["sample"] = sample
                mock["strand"] = "."
                mock.rename({"seq": "#seq"}, axis=1, inplace=True)
                mock = mock[["#seq", "start", "end", "chrom", "score", "strand", "sample"]]
                mock.sort_values(["#seq", "start", "end"], inplace=True)
                out_file = DIR_RES.joinpath(
                    "asm_label_qc", "norm_tables", "centromeres",
                    f"{sample}.active_asat_HOR_arrays_v3.bed"
                )
                mock.to_csv(out_file, sep="\t", header=True, index=False)
                mock_file = out_file.with_suffix(".bed.MOCK")
                with open(mock_file, "w"):
                    pass
    # END OF RUN BLOCK


localrules: add_centro_merge_label
rule add_centro_merge_label:
    input:
        bed = DIR_RES.joinpath(
            "asm_label_qc", "norm_tables", "centromeres",
            "{sample}.active_asat_HOR_arrays_v3.bed"
        )
    output:
        bed = DIR_RES.joinpath(
            "asm_label_qc", "merge_tables", "by-sample",
            "{sample}", "{sample}.centro.mrg-labels.bed"
        )
    run:
        import pandas as pd
        df = pd.read_csv(input.bed, sep="\t", header=0)
        if df.empty:
            with open(output.bed, "w"):
                pass
        else:
            df["length"] = (df["end"] - df["start"]).astype(int)
            df["merge_label"] = "CENTRO::" + df["length"].astype(str)
            df = df[["#seq", "start", "end", "merge_label"]]
            df.to_csv(output.bed, sep="\t", header=False, index=False)
    # END OF RUN BLOCK


rule run_all_normalize_centromeres:
    input:
        beds = rules.split_centromere_annotation.output.beds,
        mrg = expand(
            rules.add_centro_merge_label.output.bed,
            sample=SAMPLES
        )
