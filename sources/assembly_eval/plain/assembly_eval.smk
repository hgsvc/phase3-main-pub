import pandas as pd
import os
import numpy as np
import pysam
import gzip
from Bio import SeqIO
import re

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)


configfile: "config.yaml"


SAMPLES = list(config.keys())

FOFNS = {}
TYPE_MAP = {}
#GET_BED = []
SAMFLAG = config.get("SAMFLAG", "256")

expand_series = pd.Series([[], [], []], index=["sample", "tech", "type_map"])

for sample in SAMPLES:
    for tech in config[sample]["fofns"]:
        type_map = config[sample]["type_map"][tech]
        expand_series["sample"].append(sample)
        expand_series["tech"].append(tech)
        expand_series["type_map"].append(type_map)


for SM in SAMPLES:
    config[SM]["pic_per_png"] = config[SM].get("pic_per_png", 5)
    if "fofns" in config[SM]:
        FOFNS[SM] = {}
        if "HiFi" in config[SM]["fofns"]:
            FOFNS[SM]["HiFi"] = [
                line.strip() for line in open(config[SM]["fofns"]["HiFi"]).readlines()
            ]
        if "ONT" in config[SM]["fofns"]:
            FOFNS[SM]["ONT"] = [
                line.strip() for line in open(config[SM]["fofns"]["ONT"]).readlines()
            ]
    if "type_map" in config[SM]:
        TYPE_MAP[SM] = {}
        if "HiFi" in config[SM]["type_map"]:
            TYPE_MAP[SM]["HiFi"] = config[SM]["type_map"]["HiFi"]
        if "ONT" in config[SM]["type_map"]:
            TYPE_MAP[SM]["ONT"] = config[SM]["type_map"]["ONT"]
#    if config[SM]["genBed"]:
#        GET_BED.append(SM)


# print(FOFNS)
#shell.prefix("source %s/env.cfg; " % (SNAKEMAKE_DIR))


def getAsm(wildcards):
    asm_dict = {
            "h1": config[wildcards.sample]["asm_h1"],
            "h2": config[wildcards.sample]["asm_h2"],
        }
    try:
        tig_fofn = config[wildcards.sample]["extra_tigs"]
        with open(tig_fofn, 'r') as infile:
            tig_dict = {f'tig{i}' : line.rstrip() for i, line in enumerate(infile) }
        return asm_dict | tig_dict
    except KeyError:
        return asm_dict

def findReads(wildcards):
    return FOFNS[wildcards.sample][wildcards.tech][int(wildcards.read)]


def gatherAlignments(wildcards):
    SM = wildcards.sample
    tech = wildcards.tech
    type_map = TYPE_MAP[SM][tech]
    IDS = range(len(FOFNS.get(SM).get(tech)))
    bams = expand(SM + "/alignments/" + tech + "/" + type_map + "/tmp/{ID}.bam", ID=IDS)
    return bams


def getFai(wildcards):
    fai = str(config[wildcards.sample]["assembly"]) + ".fai"
    return fai


def gatherSplits(wildcards):
    SPLITS = glob_wildcards(
        "{sample}/coverage/contigs_{split}.bed".format(
            sample=wildcards.sample, split="{split}"
        )
    ).split
    return expand(
        rules.get_depth_cov_all.output.depth,
        sample=wildcards.sample,
        split=SPLITS,
        tech=wildcards.tech,
        type_map=TYPE_MAP[wildcards.sample][wildcards.tech],
    )


def getContigs(wildcards):
    return config[wildcards.sample]["regions"][wildcards.region]


def getNuc(wildcards):
    bam = f"{wildcards.sample}/alignments/{wildcards.tech}/{wildcards.type_map}/all_{wildcards.tech}.cram"
    # bam = f"{wildcards.sample}/alignments/{wildcards.tech}/{wildcards.type_map}/{wildcards.region}/all_{wildcards.tech}.bam"
    rptm = f"{wildcards.sample}/repeatMasker/{wildcards.region}.fa.out"
    if config[wildcards.sample]["repeat_mask"]:
        return {"bam": bam, "rptm": rptm}
    else:
        return {"bam": bam}


def getNucOpts(wildcards):
    if "nuc_opts" in config[wildcards.sample]:
        return config[wildcards.sample]["nuc_opts"]
    else:
        return ""


def getSpecies(wildcards):
    return config[wildcards.sample]["species"]


def gatherBreaks(wildcards):
    # print("here")
    BREAKS = glob_wildcards(
        "{sample}/nucFreq/{region}/contigs_{breaks}.bed".format(
            sample=wildcards.sample, region=wildcards.region, breaks="{breaks}"
        )
    ).breaks
    return expand(
        rules.nucFreq.output.png,
        sample=wildcards.sample,
        region=wildcards.region,
        tech=wildcards.tech,
        type_map=wildcards.type_map,
        breaks=BREAKS,
    )


def getFlag(wildcards):
    if wildcards.tech == "HiFi":
        return "map-pb"
    elif wildcards.tech == "ONT":
        return "map-ont"


def getCount(wildcards):
    return config[wildcards.sample]["pic_per_png"]


def aggregate_input(wildcards):
    file = checkpoints.filter_depth_cov.get(
        sample=wildcards.sample,
        tech=wildcards.tech,
        type_map=wildcards.type_map,
        scatteritem=wildcards.scatteritem,
    ).output[0]
    if os.stat(file).st_size == 0:
        return rules.correct_assemblies.output.correct_bed
    else:
        return rules.find_breaks.output.break_bed


def nucfeq_output(wildcards):
    nucfeq_out = []
    for sample in SAMPLES:
        if "regions" in config[sample]:
            for region in config[sample]["regions"]:
                if region != "":
                    for tech in config[sample]["fofns"]:
                        type_map = config[sample]["type_map"][tech]
                        nucfeq_out.append(
                            f"{sample}/nucFreq/{tech}/{region}/{type_map}/{sample}.done"
                        )
    return set(nucfeq_out)


scattergather:
    split=30,


localrules:
    rptm,
    nucFreq,
    calc_cov,
    trigger_steps,
    getFinal,
    all_eval,
    all_align,
    all_nucfreq,


rule all_eval:
    input:
        expand(
            "{sample}/assembly_eval/{tech}/{type_map}/final_qc.bed",
            zip,
            sample=expand_series["sample"],
            tech=expand_series["tech"],
            type_map=expand_series["type_map"],
        ),


rule all_align:
    input:
        expand(
            "{sample}/alignments/{tech}/{type_map}/all_{tech}.cram",
            zip,
            sample=expand_series["sample"],
            tech=expand_series["tech"],
            type_map=expand_series["type_map"],
        ),


rule all_nucfreq:
    input:
        nucfeq_output,


rule combine_asm:
    input:
        unpack(getAsm),
    output:
        comb_fast="{sample}/asm/both_haps.fasta.gz",
        comb_fai="{sample}/asm/both_haps.fasta.gz.fai",
    threads: 1
    resources:
        mem=8,
        load=100,
        disk=0,
        hrs=4,
    singularity:
        "docker://eichlerlab/assembly_eval:0.2"
    shell:
        """
        for asm in $( echo {input} ); do
            if ( file $( readlink -f ${{asm}} ) | grep -q zip ); then
                zcat ${{asm}} | bgzip -c >> {output.comb_fast}
            else
                cat ${{asm}} | bgzip -c >> {output.comb_fast}
            fi
        done
        samtools faidx {output.comb_fast}
        """


rule countKmers:
    input:
        assembly=rules.combine_asm.output.comb_fast,
    output:
        dirct=directory("{sample}/alignments/merylDB_{sample}"),
    threads: 1
    resources:
        mem=12,
        load=100,
        disk=0,
        hrs=8,
    singularity:
        "docker://eichlerlab/assembly_eval:0.2"
    shell:
        """
        meryl count k=15 output {output.dirct} {input.assembly}
        """


rule getRepeatKmers:
    input:
        db=rules.countKmers.output.dirct,
    output:
        rep="{sample}/alignments/repetitive_k15_{sample}.txt",
    threads: 1
    resources:
        mem=8,
        load=100,
        disk=0,
        hrs=8,
    singularity:
        "docker://eichlerlab/assembly_eval:0.2"
    shell:
        """
        meryl print greater-than distinct=0.9998 {input.db} > {output.rep}
        """


rule map_minimap:
    input:
        fasta=findReads,
        assembly=rules.combine_asm.output.comb_fast,
        fai=rules.combine_asm.output.comb_fai,
    output:
        bam=temp("{sample}/alignments/{tech}/minimap2/tmp/{read}.bam"),
    threads: 12
    params:
        winn_opt=getFlag,
    resources:
        mem=12,
        load=100,
        disk=0,
        hrs=48,
    singularity:
        "docker://eichlerlab/assembly_eval:0.2"
    shell:
        """
        minimap2 -a -t {threads} -I 10G -Y -x {params.winn_opt} --eqx -L --cs {input.assembly} {input.fasta} | samtools sort -o {output.bam} -
        """


rule map_winnowmap:
    input:
        fasta=findReads,
        assembly=rules.combine_asm.output.comb_fast,
        fai=rules.combine_asm.output.comb_fai,
        repKmers=rules.getRepeatKmers.output.rep,
    output:
        bam=temp("{sample}/alignments/{tech}/winnowmap/tmp/{read}.bam"),
    threads: 12
    resources:
        mem=lambda wildcards, attempt: attempt * 10,
        disk=0,
        load=100,
        hrs=120,
    params:
        winn_opt=getFlag,
    singularity:
        "docker://eichlerlab/assembly_eval:0.2"
    shell:
        """
        winnowmap -W {input.repKmers} -t {threads} -I 10G -Y -ax {params.winn_opt} --MD --cs -L --eqx {input.assembly} {input.fasta} | samtools sort -o {output.bam} -
        """


rule combine_alignments:
    input:
        align=gatherAlignments,
    output:
        combined=temp("{sample}/alignments/{tech}/{type_map}/all_{tech}.bam"),
    threads: 4
    resources:
        mem=8,
        disk=0,
        load=100,
        hrs=48,
    singularity:
        "docker://eichlerlab/assembly_eval:0.2"
    shell:
        """
        samtools merge -@ {threads} {output.combined} {input.align}
        samtools index {output.combined}
        """


rule cram_alignments:
    input:
        bam = rules.combine_alignments.output.combined,
        ref = rules.combine_asm.output.comb_fast,
    output:
        cram="{sample}/alignments/{tech}/{type_map}/all_{tech}.cram",
    threads: 4
    resources:
        mem=8,
        disk=0,
        load=100,
        hrs=48,
    singularity:
        "docker://eichlerlab/assembly_eval:0.2"
    shell:
        """
        samtools view -@ {threads} -C -T {input.ref} -o {output.cram} {input.bam}
        samtools index {output.cram}
        """

rule cram_nucfreq:
    input:
        bam = rules.cram_alignments.output.cram,
        ref = rules.combine_asm.output.comb_fast,
    output:
        cram="{sample}/alignments/nucfreq/{tech}/{type_map}/all_{tech}.cram",
    params:
        samflag=SAMFLAG
    threads: 4
    resources:
        mem=8,
        disk=0,
        load=100,
        hrs=48,
    singularity:
        "docker://eichlerlab/assembly_eval:0.2"
    shell:
        """
        samtools view -@ {threads} -F {params.samflag} -C -T {input.ref} -o {output.cram} {input.bam}
        samtools index {output.cram}
        """



rule hifi_fai_to_beds:
    input:
        fai=rules.combine_asm.output.comb_fai,
    output:
        bed=scatter.split("temp/{{sample}}/coverage/contigs_{scatteritem}.bed"),
    threads: 1
    resources:
        mem=12,
        disk=0,
        load=50,
        hrs=24,
    run:
        file_name = input.fai
        N_IDS = len(output.bed)
        # print(output.bed)
        outs = [open(f, "wb") for f in output.bed]
        # print(outs)
        fai_df = pd.read_csv(
            input.fai,
            header=None,
            names=["contig", "len", "offset", "byte", "encode"],
            sep="\t",
        )
        bed_df = pd.DataFrame(columns=["chrom", "start", "stop"])
        bed_df["chrom"] = fai_df["contig"].values.tolist()
        bed_df = bed_df.assign(start="0")
        bed_df["stop"] = fai_df["len"].values.tolist()
        out_idx = 0
        # print(out_idx)
        for i in bed_df.index:
            out_list = (
                bed_df.iloc[i, :].to_string(header=False, index=False).split("\n")
            )
            out_str = "".join(out_list)
            print(out_idx)
            # print(outs[out_idx])
            outs[out_idx].write((out_str + "\n").encode())
            out_idx += 1
            if out_idx == N_IDS:
                out_idx = 0
        for out in outs:
            out.close()


rule get_depth_cov:
    input:
        bed="temp/{sample}/coverage/contigs_{scatteritem}.bed",
        bam=rules.cram_nucfreq.output.cram,
    output:
        depth="{sample}/coverage/{tech}/{type_map}/tmp/filtered_{scatteritem}.tsv",
    threads: 1
    resources:
        mem=24,
        disk=0,
        load=50,
        hrs=12,
    singularity:
        "docker://eichlerlab/assembly_eval:0.2"
    shell:
        """
        samtools depth -b {input.bed} -a {input.bam} > {output.depth}
        """


rule calc_cov:
    input:
        tsv=gather.split(
            "{{sample}}/coverage/{{tech}}/{{type_map}}/tmp/filtered_{scatteritem}.tsv"
        ),
    output:
        depth="{sample}/coverage/{tech}/{type_map}/depth.txt",
    threads: 1
    resources:
        mem=80,
        load=500,
        disk=0,
        hrs=12,
    run:
        depth_list = []
        for file in input.tsv:
            depth_df = pd.read_csv(
                file, header=None, sep="\t", names=["contig", "base", "depth"]
            )
            depth_list.extend(depth_df["depth"])
        depth_avg = np.mean(depth_list)
        stand_dev = np.std(depth_list)
        total = depth_avg + stand_dev
        outfile = open(output.depth, "w+")
        outfile.write(str(total))
        outfile.close()


checkpoint filter_depth_cov:
    input:
        depth=rules.get_depth_cov.output.depth,
        cov=rules.calc_cov.output.depth,
    output:
        depth="{sample}/coverage/{tech}/{type_map}/{scatteritem}_filtered.bed",
    threads: 1
    resources:
        mem=50,
        load=50,
        disk=0,
        hrs=12,
    singularity:
        "docker://eichlerlab/assembly_eval:0.2"
    shell:
        """
        cat {input.cov}| xargs -i awk '{{if ($3 > {{}} || $3==0) printf ("%s\\t%s\\t%s\\n", $1, $2, $2)}}' {input.depth}| bedtools merge -i - > {output.depth}  
        """


rule rustybam:
    input:
        depth=rules.filter_depth_cov.output.depth,
        bam=rules.cram_nucfreq.output.cram,
    output:
        bed="{sample}/coverage/{tech}/{type_map}/{scatteritem}_intermediate.bed",
    threads: 8
    resources:
        mem=16,
        load=75,
        disk=0,
        hrs=120,
    benchmark:
        "benchmarks/rustybam_{sample}_{tech}_{type_map}_{scatteritem}.bench.txt"
    singularity:
        "docker://eichlerlab/assembly_eval:0.2"
    shell:
        """
        rustybam nucfreq  --bed {input.depth} {input.bam} > {output.bed}                
        """


rule filter_rustybam:
    input:
        bed=rules.rustybam.output.bed,
    output:
        combined="{sample}/coverage/{tech}/{type_map}/{scatteritem}_intermediate_2.bed",
    threads: 1
    resources:
        mem=50,
        load=50,
        disk=0,
        hrs=72,
    run:
        with open(input.bed, "r") as infile:
            for line in infile:
                print("Processing line:", line)
                if line.startswith("#"):
                    col_names=line.rstrip().split("\t")
                    print("Column names:", col_names)
                else:
                    break
        bed_df = pd.read_csv(input.bed, sep="\t", header=0, comment="#", names=col_names)
        #bed_df[["A", "G", "T", "C"]] = bed_df[["A", "G", "T", "C"]].apply(pd.to_numeric)
        bed_df["second_highest"] = bed_df.iloc[:,][["A", "G", "T", "C"]].apply(
            lambda row: row.nlargest(2).values[-1], axis=1
        )
        collapse_df = bed_df.loc[bed_df["second_highest"] > 5]
        collapse_df = collapse_df[["#chr", "start", "end"]]
        collapse_df = collapse_df.assign(type_err="collapse")
        print(collapse_df)
        missassembly_df = bed_df.loc[
            (bed_df["A"] == 0)
            & (bed_df["C"] == 0)
            & (bed_df["G"] == 0)
            & (bed_df["T"] == 0)
        ]
        missassembly_df = missassembly_df[["#chr", "start", "end"]]
        missassembly_df = missassembly_df.assign(type_err="misassembly")
        print(missassembly_df)
        combined_df = pd.concat([missassembly_df, collapse_df], ignore_index=True)
        combined_df = combined_df.sort_values(by=["#chr", "start", "end"])
        print(combined_df)
        combined_df.to_csv(output.combined, sep="\t", header=None, index=None)
        # misassembly_df.to_csv(output.misassembly_bed, sep="\t", header=None, index=None)



rule filter_bed:
    input:
        bed=rules.filter_rustybam.output.combined,
    output:
        filtered_bed="{sample}/coverage/{tech}/{type_map}/{scatteritem}_intermediate_3.bed",
    threads: 1
    resources:
        mem=20,
        load=50,
        disk=0,
        hrs=12,
    singularity:
        "docker://eichlerlab/assembly_eval:0.2"
    shell:
        """
        bedtools merge -i {input.bed} -c 1,4 -o count,collapse -d 5000 | awk '$4 > 2' | bedtools merge -i - -d 15000 -c 5 -o collapse | awk '{{printf "%s\\t%s\\t%s\\t%s\\n", $1, $2-5000, $3+5000,$4}}' > {output.filtered_bed}
        """


rule find_breaks:
    input:
        bed=rules.filter_bed.output.filtered_bed,
        fai=rules.combine_asm.output.comb_fai,
    output:
        break_bed="{sample}/coverage/{tech}/{type_map}/error/{scatteritem}_break.bed",
    threads: 1
    resources:
        mem=8,
        load=50,
        disk=0,
        hrs=1,
    run:
        bed_df = pd.read_csv(
            input.bed, sep="\t", header=None, names=["contig", "start", "stop", "type"]
        )
        fai_df = pd.read_csv(
            input.fai,
            sep="\t",
            header=None,
            names=["contig", "len", "offset", "byte", "encode"],
        )
        fai_subset = fai_df.loc[fai_df["contig"].isin(bed_df["contig"].tolist())]
        bed_fai = pd.merge(bed_df, fai_subset, on="contig", suffixes=["_bed", "_fai"])
        print(bed_fai)
        bed_fai_subset = bed_fai.loc[bed_fai["start"] > 25000]
        bed_fai_subset = bed_fai_subset.loc[
            bed_fai_subset["stop"] < bed_fai_subset["len"] - 25000
        ]
        print(bed_fai_subset)
        bed_fai_subset[["contig", "start", "stop", "type"]].to_csv(
            output.break_bed, sep="\t", header=None, index=None
        )


rule correct_assemblies:
    input:
        bed=rules.filter_depth_cov.output.depth,
    output:
        correct_bed=touch(
            "{sample}/coverage/{tech}/{type_map}/correct/{scatteritem}_break.bed"
        ),
    threads: 1
    resources:
        mem=4,
        load=50,
        disk=0,
        hrs=1,


rule trigger_steps:
    input:
        bed=aggregate_input,
    output:
        location_all="{sample}/coverage/{tech}/{type_map}/aggregated/{scatteritem}_break.bed",
    singularity:
        "docker://eichlerlab/assembly_eval:0.2"
    shell:
        """
        cp -l {input.bed} {output.location_all}
        """


rule combine_outputs:
    input:
        all_beds=gather.split(
            "{{sample}}/coverage/{{tech}}/{{type_map}}/aggregated/{scatteritem}_break.bed"
        ),
    output:
        cat_beds="{sample}/assembly_eval/{tech}/{type_map}/final_qc.bed",
    threads: 1
    resources:
        mem=10,
        disk=0,
        load=50,
        hrs=12,
    singularity:
        "docker://eichlerlab/assembly_eval:0.2"
    shell:
        """
        cat {input.all_beds} | awk -F'\\t' -v OFS='\\t' '{{split($4,a,","); asort(a); $4=""; delete b; for(i in a) if(!b[a[i]]++) $4=$4 a[i] ","; sub(/,$/,"",$4); print}}' > {output.cat_beds}
        """

rule formatContigs:
    input:
        contigs=getContigs,
    output:
        faidx_format=temp("{sample}/fasta/{region}.list"),
    threads: 1
    resources:
        mem=8,
        load=100,
        disk=0,
        hrs=1,
    run:
        bed_df = pd.read_csv(
            input.contigs, sep="\t", header=None, names=["contig", "start", "stop"]
        )
        outfile = open(output.faidx_format, "w+")

        for index in bed_df.index:
            contig = bed_df.at[index, "contig"]
            write_line = str(contig) + "\n"
            outfile.write(write_line)
        outfile.close()


rule filterFasta:
    input:
        assembly=rules.combine_asm.output.comb_fast,
        fai=rules.combine_asm.output.comb_fai,
        contigs=rules.formatContigs.output.faidx_format,
    output:
        filt_fasta="{sample}/fasta/{region}.fa",
    threads: 1
    resources:
        mem=12,
        load=100,
        disk=0,
        hrs=1,
    singularity:
        "docker://eichlerlab/assembly_eval:0.2"
    shell:
        """
        samtools faidx -r {input.contigs} {input.assembly} > {output.filt_fasta}
        """


rule rptm:
    input:
        assembly=rules.filterFasta.output.filt_fasta,
    output:
        repeat_out="{sample}/repeatMasker/{region}.fa.out",
    threads: 16
    log: "log/rptm/{sample}_{region}.log"
    resources:
        mem=8,
        load=1000,
        disk=0,
        hrs=12,
    params:
        directory="{sample}/repeatMasker/",
        species=getSpecies,
    singularity:
        "docker://eichlerlab/assembly_eval:0.2"
    shell:
        """
        RepeatMasker -e rmblast -species {params.species} -dir {params.directory} -pa {threads} {input.assembly} > {log} 2>&1 
        """


checkpoint splitBeds:
    input:
        contigs=getContigs,
    output:
        flag=touch("{sample}/nucFreq/{region}/tmp/.checkpoint.contigs.bed"),
        bed="{sample}/nucFreq/{region}/contigs_0.bed",
    threads: 1
    resources:
        mem=12,
        load=100,
        disk=0,
        hrs=2,
    params:
        dirct="{sample}/nucFreq/{region}/",
        img_count=getCount,
    run:
        bed_df = pd.read_csv(
            input.contigs, header=None, sep="\t", names=["contig", "start", "stop"], index_col=False
        )
        bed_df = bed_df.sample(frac=1).reset_index(drop=True)
        groups = bed_df.groupby(np.arange(len(bed_df.index)) // int(params.img_count))
        count = 0
        for splitno, split in groups:
            split.to_csv(
                params.dirct + "contigs_" + str(count) + ".bed",
                header=None,
                index=None,
                sep="\t",
            )
            count = count + 1


rule nucFreq:
    input:
        unpack(getNuc),
        bed="{sample}/nucFreq/{region}/contigs_{breaks}.bed",
    output:
        png="{sample}/nucFreq/{tech}/{region}/{type_map}/{sample}_{breaks}.output.png",
    threads: 1
    log: "log/nucFreq/{sample}_{tech}_{region}_{type_map}_{breaks}.log"
    resources:
        mem=12,
        load=100,
        disk=0,
        hrs=24,
    params:
        user_opt=getNucOpts,
    run:
        if "rptm" in input.keys():
            shell_string = (
                "{SNAKEMAKE_DIR}/scripts/NucPlot.py --bed {input.bed} -r %s {params.user_opt} %s {output.png} > {log} 2>&1"
                % (input.get("rptm"), input.get("bam"))
            )
            shell(shell_string)
        else:
            shell(
                "{SNAKEMAKE_DIR}/scripts/NucPlot.py --bed {input.bed} {params.user_opt} {input.bam} {output.png} > {log} 2>&1"
            )


rule getFinal:
    input:
        beds=gatherBreaks,
        flag=rules.splitBeds.output.flag,
    output:
        all_bed=touch("{sample}/nucFreq/{tech}/{region}/{type_map}/{sample}.done"),