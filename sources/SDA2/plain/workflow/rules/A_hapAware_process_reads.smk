# HAPLOTYPE AWARE MODE ACTIVE
rule A01H_linkAndLabelAsmHeaders:
    input:
        hap1=config["asm"],
        hap2=config["hap2"]
    output:
        hap1out= "results/A01H_hap1_asm.fasta",
        hap2out= "results/A01H_hap2_asm.fasta",
        hap1lnkU="results/A0U_hap1_asm.fasta",
        hap2lnkU="results/A0U_hap2_asm.fasta"
    params:
        workflowDir=workflow.basedir
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/A01H_linkAndLabelAsmHeaders.log"
    shell:"""
        echo "##### A01H_linkAndLabelAsmHeaders" > {log}
        echo "### Add in haplotype label to assembly headers if absent" >> {log}
        cat {input.hap1} | \
            awk 'BEGIN {{OFS="\\t"}} \
                (/^>/) \
                    {{if ($0~/haplotype1/) \
                        {{print $0}}
                    else \
                        {{print ">haplotype1-" substr($0,2)}}}} \
                (!/^>/) \
                    {{print $0}}' > {output.hap1out}
        cat {input.hap2} | \
            awk 'BEGIN {{OFS="\\t"}} \
                (/^>/) \
                    {{if ($0~/haplotype2/) \
                        {{print $0}}
                    else \
                        {{print ">haplotype2-" substr($0,2)}}}} \
                (!/^>/) \
                    {{print $0}}' > {output.hap2out}
        
        echo "### Link using universal codes" >> {log}
        ln -s {params.workflowDir}/../{output.hap1out} {params.workflowDir}/../{output.hap1lnkU} 2>> {log}
        ln -s {params.workflowDir}/../{output.hap2out} {params.workflowDir}/../{output.hap2lnkU} 2>> {log}
        """

rule A02H_faiIndexAsm:
    input:
        hap="results/A01H_hap{hapNum}_asm.fasta"
    output:
        hapfai="results/A01H_hap{hapNum}_asm.fasta.fai",
        haplnk="results/A02H_hap{hapNum}_asm.fai",
        haplnkU="results/A0U_hap{hapNum}_asm.fasta.fai"
    params:
        workflowDir=workflow.basedir
    resources:
        mem_mb=cluster_mem_mb_medium,
        cpus_per_task=cluster_cpus_per_task_baby,
        runtime=config["cluster_runtime_short"]
    conda: "../envs/sda2.main.yml"
    log: "logs/A02H_faiIndexAsm.hap{hapNum}.log"
    benchmark: "benchmark/A02H_faiIndexAsm.hap{hapNum}.tsv"
    shell:"""
    {{
        echo "##### A02H_faiIndexAsm - hap{wildcards.hapNum}" > {log}
        echo "### Create Samtools fai Index"
        samtools faidx {input.hap} --fai-idx {output.hapfai}

        echo "### Create Link" >> {log}
        ln -s {params.workflowDir}/../{output.hapfai} {params.workflowDir}/../{output.haplnk}
        echo "### Link using universal codes" >> {log}
        ln -s {params.workflowDir}/../{output.hapfai} {params.workflowDir}/../{output.haplnkU}
    }} 2>> {log}
    """

rule A10H_unionAndIndexAsms:
    input:
        hap1="results/A01H_hap1_asm.fasta",
        hap2="results/A01H_hap2_asm.fasta"
    output:
        hapX="results/A10H_asm.fasta",
        faiX="results/A10H_asm.fasta.fai",
        hapXlnkU="results/A0U_asm.fasta",
        faiXlnkU="results/A0U_asm.fasta.fai"
    params:
        workflowDir=workflow.basedir
    resources:
        mem_mb=cluster_mem_mb_medium,
        cpus_per_task=cluster_cpus_per_task_baby,
        runtime=config["cluster_runtime_short"]
    conda: "../envs/sda2.main.yml"
    log: "logs/A10H_unionAndIndexAsms.log"
    benchmark: "benchmark/A10H_unionAndIndexAsms.tsv"
    shell:"""
    {{
        echo "##### A10H_unionAndIndexAsms"
        echo "### Concat Asms" >> {log}
        cat {input.hap1} {input.hap2} > {output.hapX}

        echo "### Index Asm Union" >> {log}
        samtools faidx {output.hapX} --fai-idx {output.faiX}

        echo "### Create Link" >> {log}
        ln -s {params.workflowDir}/../{output.hapX} {params.workflowDir}/../{output.hapXlnkU}
        ln -s {params.workflowDir}/../{output.faiX} {params.workflowDir}/../{output.faiXlnkU}
    }} 2>> {log}
    """

rule A03H_mmiIndexAsm:
    input:
        hap="results/A01H_hap{hapNum}_asm.fasta"
    output:
        mmi="results/A03H_hap{hapNum}_asm.mmi"
    params:
        workflowDir=workflow.basedir,
        read_type=config["read_type"]
    resources:
        mem_mb=cluster_mem_mb_medium,
        cpus_per_task=cluster_cpus_per_task_small,
        runtime=config["cluster_runtime_short"]
    conda: "../envs/sda2.main.yml"
    log: "logs/A03H_mmiIndexAsm.hap{hapNum}.log"
    benchmark: "benchmark/A03H_mmiIndexAsm.hap{hapNum}.tsv"
    shell:"""
    {{
        echo "##### A03H_mmiIndexAsm - hap{wildcards.hapNum}" > {log}
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

# Sorts reads lexicographically by read name
rule A04H_alignReads:
    input:
        reads=lambda wildcards: readFiles[wildcards.base],
        mmi=ancient("results/A03H_hap{hapNum}_asm.mmi")
    output:
        fastq=temp("results/A04H_hap{hapNum}_aligned/A04H_{base}.fastq"),
        aligned=temp("results/A04H_hap{hapNum}_aligned/A04H_{base}.bam")
    priority: 10
    params:
        read_type=config["read_type"]
    resources:
        tmpdir=tmpDir,
        mem_mb=cluster_mem_mb_xlarge,
        cpus_per_task=cluster_cpus_per_task_large,
        runtime=config["cluster_runtime_long"]
    conda: "../envs/sda2.main.yml"
    log: "logs/A04H_alignReads.hap{hapNum}.{base}.log"
    benchmark: "benchmark/A04H_alignReads.hap{hapNum}.{base}.tsv"
    shell:"""
    {{
        echo "##### A04H_alignReads - hap{wildcards.hapNum} - bam: {wildcards.base}" > {log}
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
        tmp_collate_path=`mktemp -d -p {resources.tmpdir} A04H.collate.XXXX.tmp`
        tmp_mm2_path=`mktemp -d -p {resources.tmpdir} A04H.mm2.XXXX.tmp`
        tmp_sort_path=`mktemp -d -p {resources.tmpdir} A04H.sort.XXXX.tmp`

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
        echo "# Reads sorted by read name lexicographically" >> {log}
        minimap2 {input.mmi} {output.fastq} -a -x "$mm2_mode" -t {resources.cpus_per_task} --split-prefix "$tmp_mm2_path" | \
            samtools sort -N -T "$tmp_sort_path"/reads -m "$mem_per_cpu"M -@ "$numAdditionalThreads" -o {output.aligned}

        echo "### Delete Tmp Dirs and their Contents" >> {log}
        rm -rf "$tmp_collate_path" "$tmp_mm2_path" "$tmp_sort_path"
    }} 2>> {log}
    """

rule A05H_mergeReads:
    input:
        aln=expand("results/A04H_hap{{hapNum}}_aligned/A04H_{base}.bam", base=readFiles.keys())
    output:
        mrg=temp("results/A05H_hap{hapNum}_asm_includesAutozygousRgns.bam")
    resources:
        mem_mb=cluster_mem_mb_small,
        cpus_per_task=cluster_cpus_per_task_small,
        runtime=config["cluster_runtime_long"]
    conda: "../envs/sda2.main.yml"
    log: "logs/A05H_mergeReads.hap{hapNum}.log"
    benchmark: "benchmark/A05H_mergeReads.hap{hapNum}.tsv"
    shell:"""
    {{
        echo "##### A05H_mergeReads - hap{wildcards.hapNum}" > {log}
        numAdditionalThreads=$(echo "{resources.cpus_per_task} - 1" | bc)
        samtools merge -N {output.mrg} {input.aln} -@"$numAdditionalThreads"
    }} 2>> {log}
    """

# Rule assumes prior PG header lines exist in input bam files.
# Rule assumes lexicographic ordering of input bam files.
rule A06H_partitionReads:
    input:
        hap1_bam="results/A05H_hap1_asm_includesAutozygousRgns.bam",
        hap2_bam="results/A05H_hap2_asm_includesAutozygousRgns.bam"
    output:
        hap1_only_bam=temp("results/A06H_hap1_reads_only.bam"),
        hap2_only_bam=temp("results/A06H_hap2_reads_only.bam"),
        hap1_auto_bam=temp("results/A06H_hap1_autozygous_reads_only.bam"),
        hap2_auto_bam=temp("results/A06H_hap2_autozygous_reads_only.bam"),
        unmapped_bam="results/A06H_hapX_unmapped_reads_only.bam"
    conda: "../envs/sda2.main.yml"
    log: "logs/A06H_partitionReads.log"
    benchmark: "benchmark/A06H_partitionReads.tsv"
    params:
        workflowDir=workflow.basedir
    resources:
        mem_mb=cluster_mem_mb_small,
        cpus_per_task=cluster_cpus_per_task_small,
        runtime=config["cluster_runtime_long"]
    shell:"""
        echo "##### A06H_partitionReads" > {log}
        {params.workflowDir}/scripts/A06H_partitionReads.py {input.hap1_bam} {input.hap2_bam} {output.hap1_only_bam} {output.hap2_only_bam} {output.hap1_auto_bam} {output.hap2_auto_bam} {output.unmapped_bam} 2>> {log}
    """

rule A07H_disperseAutozygousReads:
    input:
        hap1_only="results/A06H_hap1_reads_only.bam",
        hap2_only="results/A06H_hap2_reads_only.bam",
        hap1_auto_bam="results/A06H_hap1_autozygous_reads_only.bam",
        hap2_auto_bam="results/A06H_hap2_autozygous_reads_only.bam"
    output:
        hap1_reads=temp("results/A07H_hap1_reads_unsorted.bam"),
        hap2_reads=temp("results/A07H_hap2_reads_unsorted.bam")
    conda: "../envs/sda2.main.yml"
    log: "logs/A07H_disperseAutozygousReads.log"
    benchmark: "benchmark/A07H_disperseAutozygousReads.tsv"
    params:
        workflowDir=workflow.basedir
    resources:
        mem_mb=cluster_mem_mb_small,
        cpus_per_task=cluster_cpus_per_task_small,
        runtime=config["cluster_runtime_long"]
    shell:"""
        echo "##### A07H_disperseAutozygousReads" > {log}
        {params.workflowDir}/scripts/A07H_randomlyDisperseAutozygousReads.py {input.hap1_only} {input.hap2_only} {input.hap1_auto_bam} {input.hap2_auto_bam} {output.hap1_reads} {output.hap2_reads} 2>> {log}
    """

rule A08H_sortReads:
    input:
        bam="results/A07H_hap{hapNum}_reads_unsorted.bam"
    output:
        sorted_bam=protected("results/A08H_hap{hapNum}_reads.bam"),
        sorted_bam_link="results/A0U_hap{hapNum}_reads.bam"
    params:
        read_type=config["read_type"],
        workflowDir=workflow.basedir
    priority: 10
    resources:
        tmpdir=tmpDir,
        mem_mb=cluster_mem_mb_xlarge,
        cpus_per_task=cluster_cpus_per_task_large,
        runtime=config["cluster_runtime_long"]
    conda: "../envs/sda2.main.yml"
    log: "logs/AA08H_sortReads.hap{hapNum}.log"
    benchmark: "benchmark/A08H_sortReads.hap{hapNum}.tsv"
    shell:"""
    {{
        echo "##### A08H_sortReads - hap{wildcards.hapNum}" > {log}
        echo "### Determine Node Variables" >> {log}
        mem_per_cpu="$(echo "{resources.mem_mb}/1.5/{resources.cpus_per_task}" | bc)"
        echo "Memory per cpu: $mem_per_cpu" >> {log}
        numAdditionalThreads=$(echo "{resources.cpus_per_task} - 1" | bc)

        echo "### Make Tmp Dir" >> {log}
        tmp_sort_path=`mktemp -d -p {resources.tmpdir} A08H.sort.XXXX.tmp`

        echo "### Sort Reads" >> {log}
        samtools sort -T "$tmp_sort_path"/reads -m "$mem_per_cpu"M -@ "$numAdditionalThreads" -o {output.sorted_bam} {input.bam}

        echo "### Delete Tmp Dirs and their Contents" >> {log}
        rm -rf "$tmp_sort_path"

        echo "### Link using universal codes" >> {log}
        ln -s {params.workflowDir}/../{output.sorted_bam} {params.workflowDir}/../{output.sorted_bam_link}
    }} 2>> {log}
    """

rule A09H_baiIndexBam:
    input:
        bam="results/A08H_hap{hapNum}_reads.bam"
    output:
        bai="results/A08H_hap{hapNum}_reads.bam.bai",
        lnk="results/A09H_hap{hapNum}_reads.bai",
        lnkU="results/A0U_hap{hapNum}_reads.bam.bai"
    conda: "../envs/sda2.main.yml"
    log: "logs/A09H_baiIndexBam.hap{hapNum}.log"
    benchmark: "benchmark/A09H_baiIndexBam.hap{hapNum}.tsv"
    params:
        workflowDir=workflow.basedir
    resources:
        mem_mb=cluster_mem_mb_small,
        cpus_per_task=cluster_cpus_per_task_small,
        runtime=config["cluster_runtime_long"]
    shell:"""
    {{
        echo "##### A09H_baiIndexBam - hap{wildcards.hapNum}" > {log}
        echo "### Index Bam" >> {log}
        samtools index {input.bam} -@{resources.cpus_per_task}

        echo "### Create Link" >> {log}
        ln -s {params.workflowDir}/../{input.bam}.bai {params.workflowDir}/../{output.lnk}

        echo "### Link using universal codes" >> {log}
        ln -s {params.workflowDir}/../{output.bai} {params.workflowDir}/../{output.lnkU}
    }} 2>> {log}
    """

rule temp_A:
    input:
        "results/A0U_hap1_reads.bam.bai",
        "results/A0U_hap2_reads.bam.bai"
    conda: "../envs/sda2.main.yml"
    shell: """
        # TODO Delete this rule
    """