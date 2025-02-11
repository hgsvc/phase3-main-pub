# HAPLOTYPE AWARE MODE DISABLED
rule A01_linkAsm:
    input:
        asm=config["asm"]
    output:
        asmlnk="results/A01_pri_asm.fasta",
        asmlnkU="results/A0U_pri_asm.fasta",
        asmlnkU2="results/A0U_asm.fasta"
    params:
        workflowDir=workflow.basedir
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/A01_linkAsm.log"
    shell:"""
        echo "##### A01_linkAsm" > {log}
        ln -s {input.asm} {params.workflowDir}/../{output.asmlnk} 2>> {log}

        echo "### Link using universal codes" >> {log}
        ln -s {input.asm} {params.workflowDir}/../{output.asmlnkU} 2>> {log}
        ln -s {input.asm} {params.workflowDir}/../{output.asmlnkU2} 2>> {log}
    """

rule A02_faiIndexAsm:
    input:
        asm="results/A01_pri_asm.fasta"
    output:
        fai="results/A01_pri_asm.fasta.fai",
        lnk="results/A02_asm.fai",
        lnkU="results/A0U_pri_asm.fasta.fai",
        lnkU2="results/A0U_asm.fasta.fai"
    params:
        workflowDir=workflow.basedir
    resources:
        mem_mb=cluster_mem_mb_medium,
        cpus_per_task=cluster_cpus_per_task_baby,
        runtime=config["cluster_runtime_short"]
    conda: "../envs/sda2.main.yml"
    log: "logs/A02_indexAsm.log"
    benchmark: "benchmark/A02_indexAsm.tsv"
    shell:"""
    {{
        echo "##### A02_indexAsm" > {log}
        samtools faidx {input.asm} --fai-idx {output.fai}

        echo "### Create Link" >> {log}
        ln -s {params.workflowDir}/../{output.fai} {params.workflowDir}/../{output.lnk}

        echo "### Link using universal codes" >> {log}
        ln -s {params.workflowDir}/../{output.fai} {params.workflowDir}/../{output.lnkU}
        ln -s {params.workflowDir}/../{output.fai} {params.workflowDir}/../{output.lnkU2}
    }} 2>> {log}
    """

rule A03_mmiIndexAsm:
    input:
        hap="results/A01_pri_asm.fasta"
    output:
        mmi="results/A03_asm.mmi"
    params:
        workflowDir=workflow.basedir,
        read_type=config["read_type"]
    resources:
        mem_mb=cluster_mem_mb_medium,
        cpus_per_task=cluster_cpus_per_task_small,
        runtime=config["cluster_runtime_short"]
    conda: "../envs/sda2.main.yml"
    log: "logs/A03_mmiIndexAsm.log"
    benchmark: "benchmark/A03_mmiIndexAsm.tsv"
    shell:"""
    {{
        echo "##### A03_mmiIndexAsm" > {log}
        echo "### Determine MM2 Parameters" >> {log}
        if [ {params.read_type} = "CCS" ]
        then
            mm2_mode="map-hifi"
        else
            mm2_mode="map-pb"
        fi
        echo "Minimap2 -x preset: $mm2_mode" >> {log}

        echo "### Create Minimap2 mmi Index" >> {log}
        minimap2 -x $mm2_mode -d {output.mmi} {input.hap} -t {resources.cpus_per_task}
    }} 2>> {log}
    """

rule A04_alignReads:
    input:
        reads=lambda wildcards: readFiles[wildcards.base],
        mmi=ancient("results/A03_asm.mmi")
    output:
        fastq=temp("results/A04_aligned/A04_{base}.fastq"),
        aligned=temp("results/A04_aligned/A04_{base}.bam")
    params:
        read_type=config["read_type"]
    priority: 10
    resources:
        tmpdir=tmpDir,
        mem_mb=cluster_mem_mb_xlarge,
        cpus_per_task=cluster_cpus_per_task_large,
        runtime=config["cluster_runtime_long"]
    conda: "../envs/sda2.main.yml"
    log: "logs/A04_alignReads.{base}.log"
    benchmark: "benchmark/A04_alignReads.{base}.tsv"
    shell:"""
    {{
        echo "##### A04_alignReads" > {log}
        echo "### Determine Node Variables" >> {log}
        mem_per_cpu="$(echo "{resources.mem_mb}/1.5/{resources.cpus_per_task}" | bc)"
        echo "Memory per cpu: $mem_per_cpu" >> {log}
        numAdditionalThreads=$(echo "{resources.cpus_per_task} - 1" | bc)

        echo "### Determine MM2 Parameters" >> {log}
        if [ {params.read_type} = "CCS" ]
        then
            mm2_mode="map-hifi"
        else
            mm2_mode="map-pb"
        fi
        echo "Minimap2 -x preset: $mm2_mode" >> {log}

        echo "### Make Tmp Dir" >> {log}
        tmp_collate_path=`mktemp -d -p {resources.tmpdir} A03.collate.XXXX.tmp`
        tmp_mm2_path=`mktemp -d -p {resources.tmpdir} A03.mm2.XXXX.tmp`
        tmp_sort_path=`mktemp -d -p {resources.tmpdir} A03.sort.XXXX.tmp`

        echo "### Convert to fastqs" >> {log}
        if [[ {input.reads} == *".bam" ]]
        then
            echo "### Create fastq files" >> {log}
            samtools view -h -F 2304 -u -@ "$numAdditionalThreads" {input.reads} | \
                samtools collate -O -u -@ "$numAdditionalThreads" - "$tmp_collate_path" | \
                samtools fastq -@ "$numAdditionalThreads" - > {output.fastq}
        elif [[ {input.reads} == *".fastq.gz" ]]
        then
            echo "### Unzip fastq" >> {log}
            gunzip -c {input.reads} 1> {output.fastq}
        elif [[ {input.reads} == *".fastq" ]]
        then
            cat {input.reads} 1> {output.fastq}
        elif [[ {input.reads} == *".fasta" ]]
        then
            cat {input.reads} 1> {output.fastq}
        elif [[ {input.reads} == *".fa" ]]
        then
            cat {input.reads} 1> {output.fastq}
        fi
        
        echo "### Align Reads" >> {log}
        minimap2 {input.mmi} {output.fastq} -a -x "$mm2_mode" -t {resources.cpus_per_task} --split-prefix "$tmp_mm2_path" | \
            samtools sort -T "$tmp_sort_path"/reads -m "$mem_per_cpu"M -@ "$numAdditionalThreads" -o {output.aligned}

        echo "### Delete Tmp Dirs and their Contents" >> {log}
        rm -rf "$tmp_collate_path" "$tmp_mm2_path" "$tmp_sort_path"
    }} 2>> {log}
    """

rule A05_mergeReads:
    input:
        aln=expand("results/A04_aligned/A04_{base}.bam", base=readFiles.keys())
    output:
        mrg=protected("results/A05_pri_asm_reads.bam"),
        lnkU="results/A0U_pri_reads.bam"
    params:
        workflowDir=workflow.basedir
    resources:
        mem_mb=cluster_mem_mb_small,
        cpus_per_task=cluster_cpus_per_task_small,
        runtime=config["cluster_runtime_long"]
    conda: "../envs/sda2.main.yml"
    log: "logs/A05_mergeReads.log"
    benchmark: "benchmark/A05_mergeReads.tsv"
    shell:"""
    {{
        echo "##### A05_mergeReads" > {log}
        numAdditionalThreads=$(echo "{resources.cpus_per_task} - 1" | bc)
        samtools merge {output.mrg} {input.aln} -@"$numAdditionalThreads"

        echo "### Link using universal codes" >> {log}
        ln -s {params.workflowDir}/../{output.mrg} {params.workflowDir}/../{output.lnkU}
    }} 2>> {log}
    """

rule A06_baiIndexBam:
    input:
        bam="results/A05_pri_asm_reads.bam"
    output:
        bai="results/A05_pri_asm_reads.bam.bai",
        lnk="results/A06_pri_asm_reads.bai",
        lnkU="results/A0U_pri_reads.bam.bai"
    conda: "../envs/sda2.main.yml"
    log: "logs/A06_indexBam.log"
    benchmark: "benchmark/A06_indexBam.tsv"
    params:
        workflowDir=workflow.basedir
    resources:
        mem_mb=cluster_mem_mb_small,
        cpus_per_task=cluster_cpus_per_task_small,
        runtime=config["cluster_runtime_long"]
    shell:"""
    {{
        echo "##### A06_mergeBams" > {log}
        echo "### Index Bam" >> {log}
        numAdditionalThreads=$(echo "{resources.cpus_per_task} - 1" | bc)
        samtools index {input.bam} -@"$numAdditionalThreads"

        echo "### Create Link" >> {log}
        ln -s {params.workflowDir}/../{input.bam}.bai {params.workflowDir}/../{output.lnk}

        echo "### Link using universal codes" >> {log}
        ln -s {params.workflowDir}/../{output.bai} {params.workflowDir}/../{output.lnkU}
    }} 2>> {log}
    """