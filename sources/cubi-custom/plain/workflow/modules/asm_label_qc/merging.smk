


DEFINED_MERGE_SPANS = {
    "wg": " | ",
    "ps": " | grep -v unassigned | ",
    "wg-no-ont": " | grep -v ISPCON | ",
    "ps-no-ont": " | grep -v unassigned | grep -v ISPCON | ",
}


rule merge_issue_labels:
    input:
        beds = expand(
            DIR_RES.joinpath(
                "asm_label_qc", "merge_tables", "by-sample",
                "{sample}", "{sample}.{annotation}.mrg-labels.bed"
            ),
            annotation=["sseqbreak", "flagger", "nucfreq", "merqury", "busco", "inspector", "deepvariant"],
            allow_missing=True
        )
    output:
        bed = DIR_RES.joinpath(
            "asm_label_qc", "merge_tables", "by-sample",
            "{sample}", "{sample}.merged-issues.{span}.bed.gz"
        )
    wildcard_constraints:
        span = "(" + "|".join(["wg", "ps", "wg-no-ont", "ps-no-ont"]) + ")"
    conda:
        DIR_ENVS.joinpath("bedtools.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 4096 * attempt
    params:
        grep=lambda wildcards: DEFINED_MERGE_SPANS[wildcards.span]
    shell:
        "cat {input.beds}"
        " {params.grep} "
        "sort -V -k1,1 -k2n,3n"
            " | "
        "bedtools merge -c 4 -o collapse -d -1 -i /dev/stdin"
            " | "
        "gzip > {output.bed}"


rule add_contig_size:
    input:
        table = rules.merge_issue_labels.output.bed,
        sizes = DIR_PROC.joinpath(
            "asm_label_qc", "assembly_size",
            "{sample}.sizes.txt"
        )
    output:
        table = temp(
            DIR_PROC.joinpath(
                "asm_label_qc", "augment_merge_tables", "temp",
                "{sample}", "{sample}.merged-issues.{span}.aug-lengths.bed.gz"
            )
        )
    run:
        import pandas as pd
        df = pd.read_csv(
            input.table, sep="\t", header=None,
            names=["seq", "start", "end", "labels"]
        )
        size_lut = dict()
        with open(input.sizes, "r") as table:
            for line in table:
                name, length = line.strip().split()
                size_lut[name] = int(length)

        df["seq_length"] = df["seq"].apply(lambda x: size_lut[x])
        missing = set(size_lut.keys()) - set(df["seq"].values)
        if "ps" in wildcards.span:
            missing = [seq for seq in missing if "unassigned" not in seq]
        if missing:
            add_new_rows = []
            for m in missing:
                add_new_rows.append(
                    (m, 0, size_lut[m], "no-labels", size_lut[m])
                )
            add_new_rows = pd.DataFrame.from_records(
                add_new_rows, columns=["seq", "start", "end", "labels", "seq_length"]
            )
            df = pd.concat([df, add_new_rows], axis=0, ignore_index=False)
            df.sort_values(["seq", "start", "end"], inplace=True)
        df.to_csv(output.table, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule add_telomere_flag:
    input:
        telo = expand(
            WORKDIR_EVAL.joinpath(
                "proc/70-annotate/telomeres/seqtk",
                "{sample}.asm-{asm_unit}.telo.tsv"
            ),
            asm_unit=MAIN_ASSEMBLY_UNITS,
            allow_missing=True
        ),
        table = rules.add_contig_size.output.table
    output:
        table = temp(
            DIR_PROC.joinpath(
                "asm_label_qc", "augment_merge_tables", "temp",
                "{sample}", "{sample}.merged-issues.{span}.aug-telo.bed.gz"
            )
        )
    run:
        import pandas as pd
        df = pd.read_csv(input.table, sep="\t", header=0)
        df["overlaps_telomere"] = 0
        for input_file in input.telo:
            with open(input_file, "r") as annotation:
                for line in annotation:
                    seq, start, end, length = line.strip().split()
                    select_seq = df["seq"] == seq
                    select_start = int(start) < df["end"]
                    select_end = int(end) > df["start"]
                    selector = select_seq & select_start & select_end
                    if selector.any():
                        df.loc[selector, "overlaps_telomere"] = 1
        df.to_csv(output.table, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule add_centromere_flag:
    input:
        centro = rules.add_centro_merge_label.output.bed,
        table = rules.add_telomere_flag.output.table
    output:
        table = temp(
            DIR_PROC.joinpath(
                "asm_label_qc", "augment_merge_tables", "temp",
                "{sample}", "{sample}.merged-issues.{span}.aug-centro.bed.gz"
            )
        )
    run:
        import pandas as pd
        df = pd.read_csv(input.table, sep="\t", header=0)
        df["overlaps_centromere"] = 0
        with open(input.centro, "r") as annotation:
            for line in annotation:
                if line.startswith("#"):
                    continue
                columns = line.strip().split()
                seq = columns[0]
                start = int(columns[1])
                end = int(columns[2])
                select_seq = df["seq"] == seq
                select_start = int(start) < df["end"]
                select_end = int(end) > df["start"]
                selector = select_seq & select_start & select_end
                if selector.any():
                    df.loc[selector, "overlaps_centromere"] = 1
        df.to_csv(output.table, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule merge_segdups_for_label_annotation:
    input:
        tsv = DIR_RES.joinpath(
            "asm_label_qc", "norm_tables",
            "segdups", "{sample}.sd-098.tsv.gz"
        ),
    output:
        merged = temp(
            DIR_PROC.joinpath(
                "asm_label_qc", "augment_merge_tables", "temp",
                "{sample}", "{sample}.merged-segdups.sd-098.bed.gz"
            )
        )
    conda:
        DIR_ENVS.joinpath("bedtools.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    shell:
        "zcat {input.tsv}"
            " | "
        "grep -v seq"
            " | "
        "sort -V -k1,1 -k2n,3n"
            " | "
        "bedtools merge -d -1 -i /dev/stdin"
            " | "
        "gzip > {output.merged}"


rule add_segdup_flag:
    input:
        segdups = rules.merge_segdups_for_label_annotation.output.merged,
        table = rules.add_centromere_flag.output.table
    output:
        table = temp(
            DIR_PROC.joinpath(
                "asm_label_qc", "augment_merge_tables", "temp",
                "{sample}", "{sample}.merged-issues.{span}.aug-segdup.bed.gz"
            )
        )
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt
    run:
        import pandas as pd
        import gzip
        df = pd.read_csv(input.table, sep="\t", header=0)
        df["overlaps_segdup98"] = 0
        with gzip.open(input.segdups, "rt") as annotation:
            for line in annotation:
                seq, start, end = line.strip().split()
                select_seq = df["seq"] == seq
                select_start = int(start) < df["end"]
                select_end = int(end) > df["start"]
                selector = select_seq & select_start & select_end
                if selector.any():
                    df.loc[selector, "overlaps_segdup98"] = 1
        df.to_csv(output.table, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule add_ngap_sizes:
    input:
        ngaps = WORKDIR_EVAL.joinpath(
            "results/regions", "{sample}",
            "{sample}.ngaps.bed"
        ),
        table = rules.add_segdup_flag.output.table
    output:
        table = DIR_RES.joinpath(
            "asm_label_qc", "merge_tables", "by-sample",
            "{sample}", "{sample}.merged-issues.{span}.augmented.bed.gz"
        )
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    run:
        import pandas as pd
        import collections as col
        ngaps = pd.read_csv(input.ngaps, sep="\t", header=0)
        ngaps_lut = col.Counter(
            dict((ctg, total_len) for ctg, total_len in
            ngaps.groupby("#contig")["length"].sum().items())
        )
        df = pd.read_csv(input.table, sep="\t", header=0)
        df["ngap_length"] = df["seq"].apply(lambda seq: ngaps_lut[seq])
        df.to_csv(output.table, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule run_all_merge_issues:
    input:
        bed = expand(
            rules.merge_issue_labels.output.bed,
            sample=SAMPLES,
            span=["wg", "ps", "wg-no-ont", "ps-no-ont"]
        ),
        aug = expand(
            rules.add_ngap_sizes.output.table,
            sample=SAMPLES,
            span=["wg", "ps", "wg-no-ont", "ps-no-ont"]
        )
