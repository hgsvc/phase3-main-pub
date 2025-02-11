
if REPEATMASKER_OFFLINE_SETUP:

    localrules: get_repmask_env_name
    rule get_repmask_env_name:
        output:
            env_name = DIR_PROC.joinpath(
                "70-annotate", "setup_repmask", "env_name.txt"
            )
        conda:
            DIR_ENVS.joinpath("biotools", "repmask.yaml")
        shell:
            "echo ${{CONDA_PREFIX}} > {output}"


    rule copy_repeatmasker_dfam_database:
        """This is intended for RepeatMasker <4.1.6
        with DFAM 3.7
        https://www.dfam.org/releases/Dfam_3.7/families/Dfam_curatedonly.h5.gz
        (this is the default database that is downloaded at runtime if not
        other database is present)

        Release 4.1.6+ of RepeatMasker supports a new database format
        that is split over several HDF files. Updating RepeatMasker
        for the offline setup would require changing these rules.

        NB: the name must be 'Dfam.h5', otherwise the RepeatMasker
        script won't recognize the file in the LIBDIR location.
        """
        input:
            db = DIR_GLOBAL_REF.joinpath("Dfam.h5"),
            env_name = rules.get_repmask_env_name.output.env_name
        output:
            ready = DIR_PROC.joinpath(
                "70-annotate", "setup_repmask", "repmask_setup.ready"
            )
        run:
            import pathlib as pl

            with open(input.env_name, "r") as textfile:
                env_prefix = textfile.read().strip()

            # the following should not fail
            env_prefix = pl.Path(env_prefix).resolve(strict=True)

            lib_location = env_prefix.joinpath(
                "share", "RepeatMasker", "Libraries",
                "Dfam.h5"
            )
            rsync_f2f(input.db, lib_location)

            with open(output.ready, "w") as dump:
                _ = dump.write("70-annotate::setup_repmask::copy_repeatmasker_dfam_database\n")
                _ = dump.write(get_timestamp())
        # END OF RUN BLOCK

else:

    localrules: repmask_setup_immediate_ready
    rule repmask_setup_immediate_ready:
        output:
            ready = DIR_PROC.joinpath(
                    "70-annotate", "setup_repmask", "repmask_setup.ready"
                )
        run:
            with open(output.ready, "w") as dump:
                    _ = dump.write("70-annotate::setup_repmask::repmask_setup_immediate_ready\n")
                    _ = dump.write(get_timestamp() + "\n")
        # END OF RUN BLOCK


localrules: repmask_setup_test_input
rule repmask_setup_test_input:
    input:
        ready = DIR_PROC.joinpath(
            "70-annotate", "setup_repmask", "repmask_setup.ready"
        )
    output:
        mock = DIR_PROC.joinpath(
            "70-annotate", "setup_repmask", "repmask_test_input.fasta"
        )
    run:
        import random as rand

        testseq = "".join(rand.choices(list("ACGT"), k=100))
        with open(output.mock, "w") as fasta:
            _ = fasta.write(f">testseq\n{testseq}\n")
    # END OF RUN BLOCK


rule repmask_build_database:
    input:
        fasta = rules.repmask_setup_test_input.output.mock
    output:
        rm_check = DIR_PROC.joinpath(
            "70-annotate", "setup_repmask", "repmask_run.ok"
        )
    log:
        DIR_LOG.joinpath(
            "70-annotate", "setup_repmask", "repmask_run.log"
        )
    conda:
        DIR_ENVS.joinpath("biotools", "repmask.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: 4096 * attempt,
        time_hrs = lambda wildcards, attempt: attempt
    params:
        out_dir = lambda wildcards, output: pathlib.Path(output.rm_check).with_suffix(".wd"),
        species = REPEATMASKER_SPECIES
    shell:
        "RepeatMasker -pa 2 -s -dir {params.out_dir} "
        "-species {params.species} {input.fasta} &> {log}"
            " && "
        "touch {output.rm_check}"
