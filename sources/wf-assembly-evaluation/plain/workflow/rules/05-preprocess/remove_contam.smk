
rule remove_assembly_contaminants:
    """
    TODO: get assembly units per sample (see output)
    """
    input:
        rep_adap = rules.normalize_merge_ncbi_fcs_adaptor_report.output.report,
        rep_contam = rules.normalize_merge_ncbi_fcs_contamination_report.output.report,
        mrg_fasta = rules.merge_and_tag_asm_units.output.mrg_fasta
    output:
        contam = DIR_PROC.joinpath(
            "05-preprocess", "remove_contam",
            "{sample}.contaminants.tmp.fa"
        ),
        asm_units = [
            DIR_PROC.joinpath("05-preprocess", "remove_contam", f"{{sample}}.{asm_unit}.tmp.fa")
            for asm_unit in ASSEMBLY_UNITS_NO_CONTAM
        ]
    log:
        DIR_LOG.joinpath("05-preprocess", "remove_contam", "{sample}.filter.log")
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt
    params:
        script=find_script("filter_contaminants"),
        out_pattern = lambda wildcards, output: output.contam.replace("contaminants", "asm-SEQTAG")
    shell:
        "{params.script} --input {input.mrg_fasta} --strip-tags --report "
            "--adapter-table {input.rep_adap} --contamination-table {input.rep_contam} "
            "--out-contaminants {output.contam} --out-pattern {params.out_pattern} &> {log}"


rule compress_clean_assembly_sequences:
    input:
        fasta = DIR_PROC.joinpath(
            "05-preprocess", "remove_contam",
            "{sample}.{asm_unit}.tmp.fa"
        )
    output:
        fagz = DIR_RES.joinpath(
            "assemblies", "{sample}", "{sample}.{asm_unit}.fasta.gz"
        ),
        fai = DIR_RES.joinpath(
            "assemblies", "{sample}", "{sample}.{asm_unit}.fasta.gz.fai"
        ),
        gzi = DIR_RES.joinpath(
            "assemblies", "{sample}", "{sample}.{asm_unit}.fasta.gz.gzi"
        ),
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    threads: CPU_LOW
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt,
    params:
        acc_res=lambda wildcards, output: register_result(output)
    shell:
        "rm -f {output.fagz}.EMPTY ; "
        "if [ -s {input.fasta} ] ; then "
        "{{ "
        "bgzip -c -l 6 -@ {threads} {input.fasta} > {output.fagz}"
            " && "
        "samtools faidx {output.fagz} ; "
        "}} else {{ "
        "touch {output.fagz} && touch {output.fai} "
        "&& touch {output.gzi} && touch {output.fagz}.EMPTY ; "
        "}} fi ;"


localrules: define_clean_assembly_regions
rule define_clean_assembly_regions:
    input:
        asm_tags = rules.create_sequence_tags_file.output.tagfile,
        fai_files = expand(
            rules.compress_clean_assembly_sequences.output.fai,
            asm_unit=ASSEMBLY_UNITS_NO_CONTAM,
            allow_missing=True,
        )
    output:
        tag_tig = DIR_RES.joinpath(
            "regions", "{sample}", "{sample}.uncontaminated.tag-contig.bed"
        ),
        tag_name = DIR_RES.joinpath(
            "regions", "{sample}", "{sample}.uncontaminated.tag-name.bed"
        ),
    run:
        import pandas as pd
        import pathlib as pl

        tags = []
        with open(input.asm_tags, "r") as asm_tags:
            for line in asm_tags:
                if not line.strip():
                    continue
                filename, tag = line.strip().split()
                tags.append(tag)
        assert len(tags) == len(set(tags))

        regions = []
        for fai_file in input.fai_files:
            fai_name = pl.Path(fai_file).name
            fai_tag = [tag for tag in tags if tag in fai_name]
            if len(fai_tag) != 1:
                raise ValueError(f"Cannot tag FASTA index file: {fai_name}")
            fai_tag = fai_tag[0]
            fai_regions = pd.read_csv(
                fai_file, sep="\t", header=None,
                names=["contig", "end"], usecols=[0, 1]
            )
            fai_regions["start"] = 0
            fai_regions["contig"] += f".{fai_tag}"
            regions.append(fai_regions)

        regions = pd.concat(regions, axis=0, ignore_index=False)
        regions.sort_values(["contig", "start"], inplace=True)

        with open(output.tag_tig, "w") as dump:
            _ = dump.write("#")
            regions[["contig", "start", "end"]].to_csv(
                dump, sep="\t", header=True, index=False
            )

        regions["name"] = regions["contig"].apply(lambda x: x.rsplit(".", 1)[1])
        regions["contig_drop"] = regions["contig"]
        regions["contig"] = regions["contig_drop"].apply(lambda x: x.rsplit(".", 1)[0])
        regions.drop("contig_drop", axis=1, inplace=True)

        with open(output.tag_name, "w") as dump:
            _ = dump.write("#")
            regions[["contig", "start", "end", "name"]].to_csv(
                dump, sep="\t", header=True, index=False
            )
    # END OF RUN BLOCK


# TODO fails if missing assembly unit for sample
rule run_all_clean_assembly:
    input:
        asm_units = expand(
            rules.compress_clean_assembly_sequences.output.fagz,
            sample=SAMPLES,
            asm_unit=ASSEMBLY_UNITS_NO_CONTAM
        ),
        contam = expand(
            rules.compress_clean_assembly_sequences.output.fagz,
            sample=SAMPLES,
            asm_unit=["contaminants"]
        ),
        regions_tagged = expand(
            rules.define_clean_assembly_regions.output.tag_tig,
            sample=SAMPLES
        ),
        regions_untagged = expand(
            rules.define_clean_assembly_regions.output.tag_name,
            sample=SAMPLES
        )
