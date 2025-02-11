
localrules: create_sequence_tags_file
rule create_sequence_tags_file:
    input:
        asm_seq = lambda wildcards: SAMPLE_INFOS[wildcards.sample][("asm", "all", "files")]
    output:
        tagfile = DIR_PROC.joinpath(
            "05-preprocess", "merge_tag_asm", "{sample}.asm-tags.tsv"
        )
    run:
        import io
        import pathlib as pl

        # by construction, these lists are ordered and matched
        tag_names = SAMPLE_INFOS[wildcards.sample][("asm", "all", "names")]
        seq_files = SAMPLE_INFOS[wildcards.sample][("asm", "all", "files")]
        assert len(tag_names) == len(seq_files)

        tags = io.StringIO()
        for seqfile, tag in zip(seq_files, tag_names):
            filename = pl.Path(seqfile).name
            tags.write(f"{filename}\t{tag}\n")

        with open(output.tagfile, "w") as dump:
            _ = dump.write(tags.getvalue())
    # END OF RUN BLOCK


rule merge_and_tag_asm_units:
    """The merge and tag code path
    largely exists to realize the
    NCBI FCS screening in a more
    efficient manner. Combining
    all assembly units and screening
    them all in a single run implies
    that the 400 GB database only
    has to be loaded once per assembly.
    """
    input:
        asm_seq = lambda wildcards: SAMPLE_INFOS[wildcards.sample][("asm", "all", "files")],
        tags = DIR_PROC.joinpath(
            "05-preprocess", "merge_tag_asm", "{sample}.asm-tags.tsv"
        ),
        skip = lambda wildcards: SAMPLE_INFOS[wildcards.sample].get("skip_seqs", [])
    output:
        mrg_fasta = DIR_PROC.joinpath(
            "05-preprocess", "merge_tag_asm", "{sample}.asm-mrg-tag.fasta"
        )
    log:
        DIR_LOG.joinpath(
            "05-preprocess", "merge_tag_asm", "{sample}.asm-mrg-tag.log"
        )
    conda:
        DIR_ENVS.joinpath("pyutils.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt,
        time_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        script = find_script("fasta_tag_merge"),
        buffer = int(5e8),
        set_skip_arg = lambda wildcards, input: (
            f"--skip {input.skip}"
            if not isinstance(input.skip, list) and pathlib.Path(input.skip).is_file()
            else ""
        )
    shell:
        "{params.script} --input {input.asm_seq} --seq-tags {input.tags} {params.set_skip_arg} "
            "--report --buffer-size {params.buffer} --output {output.mrg_fasta} 2> {log}"


rule index_merged_tagged_assembly_fasta:
    input:
        mrg_fasta = rules.merge_and_tag_asm_units.output.mrg_fasta
    output:
        fai = DIR_PROC.joinpath(
            "05-preprocess", "merge_tag_asm", "{sample}.asm-mrg-tag.fasta.fai"
        )
    conda:
        DIR_ENVS.joinpath("biotools", "utils.yaml")
    shell:
        "samtools faidx {input.mrg_fasta}"


localrules: dump_clean_assembly_regions
rule dump_clean_assembly_regions:
    """This rule is the counterpart of
    'define_clean_assembly_regions' in the
    'remove_contam.smk' module to be executed
    in case the initial contamination scanning
    is skipped.
    """
    input:
        fai = rules.index_merged_tagged_assembly_fasta.output.fai,
        skip = lambda wildcards: SAMPLE_INFOS[wildcards.sample].get("skip_seqs", "")
    output:
        bed = DIR_PROC.joinpath(
            "05-preprocess", "merge_tag_asm", "{sample}.sequences.bed"
        )
    run:
        import pathlib as pl
        import pandas as pd
        if pl.Path(input.skip).is_file():
            with open(input.skip) as listing:
                skip_seqs = set(listing.read().strip().split())
        else:
            skip_seqs = set()
        df = pd.read_csv(
            input.fai, sep="\t", usecols=[0,1],
            header=None, names=["contig", "end"]
        )
        df["start"] = 0
        df = df.loc[~df["contig"].isin(skip_seqs), :].copy()
        df = df[["contig", "start", "end"]]
        df.sort_values("contig", inplace=True)
        with open(output.bed, "w") as dump:
            _ = dump.write("#")
            df.to_csv(dump, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


localrules: create_sed_replacement_files
rule create_sed_replacement_files:
    input:
        fai = rules.index_merged_tagged_assembly_fasta.output.fai
    output:
        tag_to_untag = DIR_RES.joinpath(
            "auxiliary", "contig_renaming",
            "{sample}.mrg-asm.tagged-to-untagged.sed"
        ),
        untag_to_tag = DIR_RES.joinpath(
            "auxiliary", "contig_renaming",
            "{sample}.mrg-asm.untagged-to-tagged.sed"
        )
    run:
        import contextlib as ctl

        with ctl.ExitStack() as exs:
            fai_file = exs.enter_context(open(input.fai, "r"))
            out_tag_to_untag = exs.enter_context(open(output.tag_to_untag, "w"))
            out_untag_to_tag = exs.enter_context(open(output.untag_to_tag, "w"))

            for line in fai_file:
                tagged_seq = line.strip().split()[0]
                untagged_seq = tagged_seq.rsplit(".", 1)[0]
                _ = out_tag_to_untag.write(
                    f"s/{tagged_seq}/{untagged_seq}/g\n"
                )
                _ = out_untag_to_tag.write(
                    f"s/{untagged_seq}/{tagged_seq}/g\n"
                )
    # END OF RUN BLOCK


rule run_merge_tag_all_assemblies:
    input:
        mrg_fasta = expand(
            rules.merge_and_tag_asm_units.output.mrg_fasta,
            sample=SAMPLES
        )
