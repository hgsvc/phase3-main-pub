localrules: load_verkko_testdata
rule load_verkko_testdata:
    output:
        hifi = DIR_LOCAL_REF.joinpath("hifi.fastq.gz"),
        nano = DIR_LOCAL_REF.joinpath("ont.fastq.gz"),
    shell:
        "curl -L https://obj.umiacs.umd.edu/sergek/shared/ecoli_hifi_subset24x.fastq.gz -o {output.hifi}"
            " && "
        "curl -L https://obj.umiacs.umd.edu/sergek/shared/ecoli_ont_subset50x.fastq.gz -o {output.nano}"


rule run_verkko_test_local:
    input:
        hifi = rules.load_verkko_testdata.output.hifi,
        nano = rules.load_verkko_testdata.output.nano,
    output:
        done = DIR_PROC.joinpath("testdata/verkko/local/assembly.ok")
    log:
        DIR_LOG.joinpath("testdata/verkko/local/assembly.log")
    benchmark:
        DIR_RSRC.joinpath("testdata/verkko/local/assembly.log")
    conda:
        DIR_ENVS.joinpath("verkko.yaml")
    threads: 6
    resources:
        mem_mb = lambda wildcards, attempt: 16384 * attempt,
        mem_gb = lambda wildcards, attempt: 16 * attempt,
        time_hrs = lambda wildcards, attempt: attempt * attempt,
    params:
        dryrun = "--dry-run" if VERKKO_DRY_RUN else "",
        check = lambda wildcards, output: "" if VERKKO_DRY_RUN else f" && touch {output.done}",
        wd=lambda wildcards, output: pathlib.Path(output.done).with_suffix(".wd")
    shell:
        "/usr/bin/time -v "
        "verkko --local "
            "--local-memory {resources.mem_gb} "
            "--local-cpus {threads} "
            "--threads {threads} "
            "--python `which python` "
            "--mbg `which MBG` "
            "--graphaligner `which GraphAligner` "
            "--mashmap `which mashmap` "
            "--hifi {input.hifi} "
            "--nano {input.nano} "
            "-d {params.wd} "
            "--snakeopts \"--restart-times 1 {params.dryrun}\" "
            " &> {log} {params.check}"


localrules: run_verkko_test_cluster
rule run_verkko_test_cluster:
    """
    This rule sets the Verkko option
    --lsf
    to force-run Verkko in some cluster
    mode. It is vital for proper execution
    that ALL Snakemake profile parameters
    are properly set in the profile supplied
    via --snakeopts to override the defaults.
    """
    input:
        hifi = rules.load_verkko_testdata.output.hifi,
        nano = rules.load_verkko_testdata.output.nano,
        profile = ancient(config["verkko_smk_profile"]),
    output:
        done = DIR_PROC.joinpath("testdata/verkko/cluster/assembly.ok")
    log:
        DIR_LOG.joinpath("testdata/verkko/cluster/assembly.log")
    benchmark:
        DIR_RSRC.joinpath("testdata/verkko/cluster/assembly.log")
    conda:
        DIR_ENVS.joinpath("verkko.yaml")
    params:
        dryrun = "--dry-run" if VERKKO_DRY_RUN else "",
        check = lambda wildcards, output: "" if VERKKO_DRY_RUN else f" && touch {output.done}",
        wd=lambda wildcards, output: pathlib.Path(output.done).with_suffix(".wd")
    shell:
        "/usr/bin/time -v "
        "verkko --lsf "
            "--python `which python` "
            "--mbg `which MBG` "
            "--graphaligner `which GraphAligner` "
            "--mashmap `which mashmap` "
            "--hifi {input.hifi} "
            "--nano {input.nano} "
            "-d {params.wd} "
            "--snakeopts \"--profile $PWD/{input.profile} {params.dryrun}\" "
            "&> {log} {params.check}"


localrules: run_verkko_tests
rule run_verkko_tests:
    input:
        assm_dirs = [
            DIR_PROC.joinpath("testdata/verkko/local/assembly.ok"),
            DIR_PROC.joinpath("testdata/verkko/cluster/assembly.ok")
        ]


rule extract_meryl_hapmer_db:
    """Rule exist to support gzipped meryl
    k-mer / hap-mer databases as input
    """
    input:
        targz = lambda wildcards: MAP_SAMPLE_TO_INPUT_FILES[wildcards.sample][f"{wildcards.hap_db}_tar"]
    output:
        meryl = directory(DIR_PROC.joinpath(
            "00-prepare", "verkko", "hapmer_dbs",
            "{sample}.{hap_db}.meryl"
        ))
    params:
        folder = lambda wildcards, output: pathlib.Path(output.meryl),
        acc_in=lambda wildcards, input: register_input(input),
    shell:
        "mkdir -p {output.meryl}"
            " && "
        "tar xzf {input.targz} --strip-components=1 -C {output.meryl}"
