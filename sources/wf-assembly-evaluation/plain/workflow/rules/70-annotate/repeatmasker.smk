
rule create_plain_assembly_file:
    """ 2023-12-19
    This rule was added to switch to the Dfam TE Tools
    container for RepeatMasker and HMMER/nhmmer (see module).
    as the default due to the possibility of a stable
    offline deployment (default Dfam database included).
    https://github.com/Dfam-consortium/TETools

    This rule simply decouples unzipping the fasta file
    and can thus make use of Snakemake's temp() wrapper.
    """
    input:
        fasta = get_asm_unit
    output:
        tmp_fa = temp(
            str(DIR_PROC.joinpath(
                "70-annotate", "repeatmasker",
                "{sample}.{asm_unit}.repmask.tmp.fa"
            ))
        )
    resources:
        time_hrs=lambda wildcards, attempt: attempt * attempt
    shell:
        "gzip -dc {input.fasta} > {output.tmp_fa}"


rule repeatmasker_assembly_run:
    """
    Uses default RepeatMasker library (Dfam/curated)
    and is designed for a Conda install.
    The TE tool container has caused problems of unclear origin.
    cf.: github.com/rmhubley/RepeatMasker/issues/233

    NB: RepeatMasker cannot process compressed files, but the
    main process does not abort, just the I/O part seems to die.
    Hence the unzip rule above.
    """
    input:
        rm_db_ok = rules.repmask_build_database.output.rm_check,
        fasta = rules.create_plain_assembly_file.output.tmp_fa
    output:
        repmask_out = multiext(
            str(DIR_PROC.joinpath(
                "70-annotate", "repeatmasker",
                "{sample}.{asm_unit}.repmask.wd",
                "{sample}.{asm_unit}.repmask.tmp.fa"
            )),
            ".masked", ".out", ".tbl"
        )
    log:
        DIR_LOG.joinpath(
            "70-annotate", "repeatmasker",
            "{sample}.{asm_unit}.repmask.log"
        )
    benchmark:
        DIR_RSRC.joinpath(
            "70-annotate", "repeatmasker",
            "{sample}.{asm_unit}.repmask.rsrc"
        )
    conda:
        DIR_ENVS.joinpath("biotools", "repmask.yaml")
    threads: CPU_LOW
    resources:
        mem_mb = lambda wildcards, attempt, input: attempt * get_repeatmasker_run_memory_mb(input.size_mb, compressed=True),
        time_hrs = lambda wildcards, attempt, input: attempt * get_repeatmasker_run_time_hrs(input.size_mb, compressed=True),
    params:
        out_dir = lambda wildcards, output: pathlib.Path(output.repmask_out[0]).parent,
        species = REPEATMASKER_SPECIES
    shell:
        "RepeatMasker -pa {threads} -s -dir {params.out_dir} "
        "-species {params.species} {input.fasta} &> {log}"


rule collect_repeatmasker_output:
    """DEBUG / FIX 2023-11-29
    It seems like RepeatMasker is only gzipping the ".cat"
    output file if it is larger than some threshold, or
    otherwise failing w/o proper signalling. In any case,
    the workflow cannot complete for a subset of runs if
    the file ".cat.gz" is an expected output. Hence, this
    rule now compensates for that total headache and searches
    for ".cat" or ".cat.gz" in the output folder of the
    RepeatMasker run and adds that file to the tar output.
    """
    input:
        repmask_files = rules.repeatmasker_assembly_run.output.repmask_out
    output:
        tar = DIR_RES.joinpath(
                "annotations", "repeatmasker",
                "{sample}.{asm_unit}.repmask.raw.tar.gz",
            )
    resources:
        time_hrs=lambda wildcards, attempt: max(0, attempt - 1)
    run:
        import pathlib as pl
        import subprocess as sp

        tar_cmd = f"tar czf {output.tar} -C "
        for idx, filepath in enumerate(input.repmask_files, start=0):
            if idx == 0:
                tar_dir = pl.Path(filepath).parent
                tar_cmd += f"{tar_dir} "
            filename = pl.Path(filepath).name
            tar_cmd += f"./{filename} "
        # see docstring above for the following stunt ...
        cat_output = list(tar_dir.glob("*.cat*"))
        if len(cat_output) != 1:
            raise RuntimeError(
                f"Incomplete RepeatMasker run for {tar_dir} - '.cat*' output missing."
            )
        else:
            cat_output = cat_output[0]
            cat_output_filename = cat_output.name
            tar_cmd += f"./{cat_output_filename}"
        try:
            _ = sp.check_output(tar_cmd, shell=True)
        except sp.CalledProcessError as spe:
            logerr(
                (
                    "tar of RepeatMasker output failed: "
                    f"{wildcards.sample} / {wildcards.asm_unit} "
                    "\nOriginal command string:\n"
                    f"{tar_cmd}"
                ))
            raise spe
    # END OF RUN BLOCK


rule normalize_repeatmasker_table:
    input:
        table = DIR_PROC.joinpath(
            "70-annotate", "repeatmasker",
            "{sample}.{asm_unit}.repmask.wd",
            "{sample}.{asm_unit}.repmask.tmp.fa.out"
        ),
    output:
        tsv = DIR_RES.joinpath(
            "annotations", "repeatmasker",
            "{sample}.{asm_unit}.repmask.matches.tsv.gz"
        ),
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt
    run:
        import gzip
        import pandas as pd

        names_explained = [
            (0, 'sw_score', 'Smith-Waterman score of the match'),
            (1, 'pct_subs', 'Pct. substitutions in matching region compared to the consensus'),
            (2, 'pct_del', 'Pct. of bases opposite a gap in the query sequence (deleted bp)'),
            (3, 'pct_ins', 'Pct. of bases opposite a gap in the repeat consensus (inserted bp)'),
            (4, 'query_name', 'Query sequence name'),
            (5, 'query_start', 'Starting position of match in query sequence'),
            (6, 'query_end', 'Ending position of match in query sequence'),
            (7, 'query_after_bp', 'bp in query sequence past the ending position of match'),
            (8, 'is_complement_match', 'True: match is with the Complement of the consensus sequence in the database'),
            (9, 'repeat_name', 'Name of the matching interspersed repeat'),
            (10, 'repeat_class', 'Class of the repeat'),
            (11, 'repeat_start', 'Starting position of match in database sequence (using top-strand numbering)'),
            (12, 'repeat_end', 'Ending position of match in database sequence'),
            (13, 'repeat_before_bp', 'bp in (complement of) the repeat consensus sequence prior to beginning of the match'),
            (14, 'match_ID', 'Consecutive ID of match'),
            (15, 'is_partly_included', 'True: there is a higher-scoring match whose domain partly (<80%) includes the domain of this match')
        ]

        df = pd.read_csv(
            input.table,
            header=None,
            index_col=False,
            delimiter=r"\s+",
            skip_blank_lines=True,
            skiprows=[0,1],
            comment='#',
            names=[n[1] for n in names_explained]
        )

        for c in ['query_after_bp', 'repeat_start', 'repeat_before_bp']:
            df[c] = df[c].apply(lambda x: int(x.strip('()')))
            df[c] = df[c].astype(int)

        df['is_complement_match'] = df['is_complement_match'].apply(lambda x: True if x == 'C' else False)
        df['is_complement_match'] = df['is_complement_match'].astype(bool)
        df['is_partly_included'] = df['is_partly_included'].apply(lambda x: True if x == '*' else False)
        df['is_partly_included'] = df['is_partly_included'].astype(bool)
        df.sort_values(['query_name', 'query_start', 'sw_score'], ascending=[True, True, False], inplace=True)

        with gzip.open(output.tsv, 'wt') as table:
            for pos, head, comment in names_explained:
                _ = table.write(f'##.{pos} - {head} - {comment}\n')
        df.to_csv(output.tsv, mode="a", sep='\t', header=True, index=False)
    # END OF RUN BLOCK


rule run_all_repeatmasker_jobs:
    """NB: failed RepeatMasker runs do not clean
    up after themselves (presumably to keep debugging
    information intact) and there is no switch to change
    that behavior. Hence, this trigger rule also performs
    the cleanup operation. Obviously, this can only be
    kicked off after all jobs have passed (restarted until
    completion), which makes this a typical candidate for
    leaving behind garbage in case the pipeline is
    interrupted in some way. Extremely annoying!!!
    """
    input:
        norm_tables = expand(
            rules.normalize_repeatmasker_table.output.tsv,
            sample=SAMPLES,
            asm_unit=ASSEMBLY_UNITS_NO_CONTAM
        ),
        tars = expand(
            rules.collect_repeatmasker_output.output.tar,
            sample=SAMPLES,
            asm_unit=ASSEMBLY_UNITS_NO_CONTAM
        ),
    shell:
        "rm -rf RM_*"

