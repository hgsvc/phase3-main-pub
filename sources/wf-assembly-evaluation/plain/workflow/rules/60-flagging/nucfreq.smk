
rule build_nucfreq_cache:
    input:
        script = rules.clone_nucfreqtwo_repo.output.nucfreqtwo,
        bam = rules.merge_read_to_assembly_subset_alignments.output.bam,
        bai = rules.merge_read_to_assembly_subset_alignments.output.bai
    output:
        hdf = DIR_PROC.joinpath(
            "60-flagging", "nucfreq",
            "{sample}.{read_type}.{aln_subset}.nfdata.hdf"
        )
    benchmark:
        DIR_RSRC.joinpath(
            "60-flagging", "nucfreq",
            "{sample}.{read_type}.{aln_subset}.nucfreqtwo.rsrc"
        )
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    threads: CPU_LOW
    resources:
        time_hrs = lambda wildcards, attempt: 8 * attempt,
        mem_mb = lambda wildcards, attempt: 24576 + 24576 * attempt
    shell:
        "{input.script} --infile {input.bam} --out-hdf-cache {output.hdf} "
        "--mode process --threads {threads} "
        "--flag-discordant-pct 10 "
        "--flag-discordant-abs 2 "
        "--flag-min-interval 500 "
        "--flag-num-hets 5 "
        "--flag-store-window 50000 "
        "--min-region-size 500000"


rule dump_nucfreq_flagged_regions:
    input:
        hdf = expand(
            rules.build_nucfreq_cache.output.hdf,
            read_type="hifi",
            aln_subset="onlyPRI",
            allow_missing=True
        )
    output:
        bed = DIR_RES.joinpath(
            "regions", "{sample}",
            "{sample}.nucfreq-flagged.bed.gz"
        )
    run:
        import pandas as pd
        import gzip

        column_sort = [
            "contig", "start", "end",
            "asm_unit", "num_hets", "het_pct"
        ]
        # due to expand, turned into a list
        hdf_cache = input.hdf[0]

        dump_regions = []
        with pd.HDFStore(hdf_cache, "r") as hdf:
            contig_map = hdf["/group_table"]
            for key in hdf.keys():
                if "regions" not in key:
                    continue
                this_group = int(key.split("/")[1].strip("group_"))
                this_contig = contig_map.loc[contig_map["group"] == this_group, "contig"]
                assert this_contig.shape[0] == 1
                this_contig, asm_unit = this_contig.values[0].rsplit(".", 1)

                data = hdf[key]
                data["contig"] = this_contig
                data["asm_unit"] = asm_unit
                data.reset_index(drop=True, inplace=True)
                dump_regions.append(data)

        dump_regions = pd.concat(dump_regions, axis=0, ignore_index=False)
        dump_regions.sort_values(["contig", "start", "end"], inplace=True)
        dump_regions = dump_regions[column_sort]

        with gzip.open(output.bed, "wt") as dump:
            _ = dump.write("#")
            dump_regions.to_csv(dump, header=True, index=False, sep="\t")
    # END OF RUN BLOCK


rule annotate_nucfreq_regions_readcov:
    input:
        regions = rules.dump_nucfreq_flagged_regions.output.bed,
        readcov = rules.transform_mosdepth_window_read_coverage.output.hdf
    output:
        tsv = DIR_RES.joinpath(
            "regions", "{sample}",
            "{sample}.nucfreq.covann.tsv.gz"
        )
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    params:
        script=find_script("annotate_region_cov")
    shell:
        "{params.script} --regions {input.regions} --read-cov {input.readcov} "
        "--output {output.tsv}"


rule run_all_nucfreq_jobs:
    input:
        hdf = expand(
            rules.build_nucfreq_cache.output.hdf,
            sample=SAMPLES,
            read_type=["hifi"],
            aln_subset=["onlyPRI"]
        ),
        regions = expand(
            rules.annotate_nucfreq_regions_readcov.output.tsv,
            sample=SAMPLES
        ),
