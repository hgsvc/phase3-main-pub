
rule run_ribotin_on_verkko:
    input:
        verkko_asm_check = WORKDIR_ASSEMBLY.joinpath(
            "proc/10-assemble/verkko", "{sample}.ps-none.ok"
        )
    output:
        checkfile = DIR_PROC.joinpath(
            "rdna", "ribotin", "{sample}.ps-none.ribotin.ok"
        )
    benchmark:
        DIR_RSRC.joinpath(
            "rdna", "ribotin", "{sample}.ps-none.ribotin.ok"
        )
    conda:
        DIR_ENVS.joinpath("ribotin.yaml")
    threads: 6
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 6144,
        time_hrs=lambda wildcards, attempt: attempt * 11
    params:
        asm_folder = lambda wildcards, input: pathlib.Path(input.verkko_asm_check).with_suffix(".wd"),
        out_folder = lambda wildcards, output: pathlib.Path(output.checkfile).with_suffix(".wd"),
        tmp_folder = lambda wildcards, output: pathlib.Path(output.checkfile).with_suffix(".wd").joinpath("tmp"),
    shell:
        "ribotin-verkko --sample-name {wildcards.sample} -x human "
        "-i {params.asm_folder} -o {params.out_folder} "
        "--mbg `which MBG` --graphaligner `which GraphAligner` "
        "--ul-tmp-folder {params.tmp_folder} "
            " && "
        "touch {output.checkfile}"


## DEBUG
# sample NA19331 fails - diagnosed by Mikko as follows:
# (via Slack / Oct. 26. 2023)
#
# Looks like the genome has some weird repeats going on within the rDNA which break ribotin's assumptions.
# There's a repeat type that the KY962518.1 rDNA reference calls "similar to Long Repeat 1" which usually
# seems to be a few kbp long, but in this sample on one (or more) of the morphs it's at least 30kbp long
# so too long to be spanned by hifi reads. This means the resulting MBG graph is no longer locally acyclic which ribotin assumed.
# There's no quick fix so for now I'd say just leave this sample out of the analysis until I've fixed this in a future version

DEBUG_RIBOTIN_SAMPLES = [smp for smp in PLAIN_SAMPLES if smp != "NA19331"]

rule run_all_ribotin:
    input:
        checkfiles = expand(
            rules.run_ribotin_on_verkko.output.checkfile,
            sample=DEBUG_RIBOTIN_SAMPLES
        )
