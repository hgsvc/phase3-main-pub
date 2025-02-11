# Flow of this smk file:
# - D01 map asm to model (getting resolved originals)
# - D02 filter keeping >50% aligned resolved original hits to seq (and index bam) (out=bam)
# - D03 bamToBed (out bed)
# - D04 filter out single exon genes
# - D05 filter >= min_hit_length
# - D06 get fasta of resolved original hits
# - D07 properly name resolved original hits fasta


rule D01_FindResolvedOriginals:
    input:
        gm="results/C03_gene_model_filt.fasta",
        asm="results/A0U_asm.fasta"
    output:
        bam="results/D01_resolved_originals_prefilt.bam",
        csi="results/D01_resolved_originals_prefilt.bam.csi"
    resources:
        mem_mb=cluster_mem_mb_xlarge,
        cpus_per_task=cluster_cpus_per_task_medium,
        runtime=config["cluster_runtime_long"],
        tmpdir=tmpDir
    conda: "../envs/sda2.main.yml"
    log: "logs/D01_FindResolvedOriginals.log"
    benchmark: "benchmark/D01_FindResolvedOriginals.tsv"
    shell:"""
    {{
        echo "##### D01_FindResolvedOriginals" > {log}
        echo "### Determine Node Variables" >> {log}
        mem_per_cpu="$(echo "{resources.mem_mb}/1.5/{resources.cpus_per_task}" | bc)"
        echo "Memory per cpu: $mem_per_cpu" >> {log}

        echo "### Make Tmp Dir" >> {log}
        tmp_sort_path=`mktemp -d -p {resources.tmpdir} geneModel.sort.XXXX.tmp`

        echo "### Minimap Gene Model" >> {log}
        minimap2 -x splice -a -t {resources.cpus_per_task} {input.asm} {input.gm} | \
            samtools sort -T "$tmp_sort_path"/originalHits -m "$mem_per_cpu"M -@ $(( {resources.cpus_per_task}-1 )) -o {output.bam}

        echo "### Delete Samtools Sort Tmp Dir and its Contents" >> {log}
        tmp_file_size="$( du -sh $tmp_sort_path )"
        echo "Samtools sort temp directory size: $tmp_file_size ($tmp_sort_path)" >> {log}
        rm -rf "$tmp_sort_path"

        echo "### Index Resolved Original Hits Bam" >> {log}
        samtools index -c {output.bam}
    }} 2>> {log}
    """

rule D02_FilterResolvedOriginals_FiltPercentAligned:
    input:
        bam="results/D01_resolved_originals_prefilt.bam"
    output:
        filt="results/D02_resolved_originals_filtPercentAligned.bam",
        idx="results/D02_resolved_originals_filtPercentAligned.bam.csi"
    params:
        min_gm_alignment=config["min_gene_model_alignment"],
        workflowDir=workflow.basedir
    resources:
        mem_mb=cluster_mem_mb_baby,
        cpus_per_task=cluster_cpus_per_task_baby,
        runtime=config["cluster_runtime_short"]
    conda: "../envs/sda2.main.yml"
    log: "logs/D02_FilterResolvedOriginals_FiltPercentAligned.log"
    shell:"""
    {{
        echo "##### D02_FilterResolvedOriginals_FiltPercentAligned" > {log}
        echo "### Filter Gene Model" >> {log}
        {params.workflowDir}/scripts/D02_AnnotateMappedLength.py {input.bam} | \
            samtools view -e "[pa] >= {params.min_gm_alignment}" -b -o {output.filt}

        echo "### Index Resolved Original Hits Filt Bam" >> {log}
        samtools index -c {output.filt}
    }} 2>> {log}
    """

rule D03_BamToBed:
    input:
        bam="results/D02_resolved_originals_filtPercentAligned.bam"
    output:
        bed="results/D03_resolved_originals_filtPercentAligned.bed"
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/D03_BamToBed.log"
    benchmark: "benchmark/D03_BamToBed.tsv"
    shell:"""
        echo "##### D03_BamToBed" > {log}
        bedtools bamtobed -bed12 -i {input.bam} 1> {output.bed} 2>> {log}
    """

rule D04_FilterResolvedOriginals_FiltMultiExon:
    input:
        bed="results/D03_resolved_originals_filtPercentAligned.bed"
    output:
        filt="results/D04_resolved_originals_filtMultiExon.bed"
    params:
        filt_out_single_exons=config["flag_filt_single_exon_genes"],
        workflowDir=workflow.basedir
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/D04_FilterResolvedOriginals_FiltMultiExon.log"
    shell:"""
    {{
        echo "##### D04_FilterResolvedOriginals_FiltMultiExon" > {log}
        if [ {params.filt_out_single_exons} = "True" ]
        then
            echo "### Remove Genes with only a single exon" >> {log}
            cat {input.bed} | \
                awk ' \
                    BEGIN \
                        {{OFS="\\t"}} \
                    ($10 > 1) \
                        {{print}}' 1> {output.filt}
        else
            echo "### Do not remove genes with only a single exon: create symlink instead." >> {log}
            ln -s {params.workflowDir}/../{input.bed} {params.workflowDir}/../{output.filt}
        fi
    }} 2>> {log}
    """

rule D05_FilterResolvedOriginals_FiltMinLength:
    input:
        bed="results/D04_resolved_originals_filtMultiExon.bed"
    output:
        filt="results/D05_resolved_originals_filtMinLength.bed"
    params:
        min_hit_length=config["min_hit_length"]
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/D05_FilterResolvedOriginals_FiltMinLength.log"
    benchmark: "benchmark/D05_FilterResolvedOriginals_FiltMinLength.tsv"
    shell:"""
        echo "##### D05_FilterResolvedOriginals_FiltMinLength" > {log}
        cat {input.bed} | \
            awk -v minLength={params.min_hit_length} ' \
                BEGIN \
                    {{OFS="\\t"}} \
                ($3-$2 >= minLength) \
                    {{print}}' 1> {output.filt} 2>> {log}
    """

rule D06_GetResolvedOriginalsFasta_unnamed:
    input:
        bed="results/D05_resolved_originals_filtMinLength.bed",
        asm="results/A0U_asm.fasta"
    output:
        rgn=temp("results/D06_resolved_originals_filtMinLength.rgn"),
        fa="results/D06_resolved_originals_unnamed.fasta"
    resources:
        mem_mb=cluster_mem_mb_baby,
        cpus_per_task=cluster_cpus_per_task_baby,
        runtime=config["cluster_runtime_short"]
    conda: "../envs/sda2.main.yml"
    log: "logs/D06_GetResolvedOriginalsFasta_unnamed.log"
    benchmark: "benchmark/D06_GetResolvedOriginalsFasta_unnamed.tsv"
    shell:"""
        echo "##### D06_GetResolvedOriginalsFasta_unnamed" > {log}
        cat {input.bed} | awk '{{print $1":"$2"-"$3 }}' 1> {output.rgn} 2>> {log}
        samtools faidx {input.asm} -r {output.rgn} 1> {output.fa} 2>> {log}
    """

rule D07_GetResolvedOriginalsFasta:
    input:
        fa_unnamed="results/D06_resolved_originals_unnamed.fasta",
        bed="results/D05_resolved_originals_filtMinLength.bed"
    output:
        fa="results/D07_resolved_originals.fasta"
    params:
        workflowDir=workflow.basedir
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/D07_GetResolvedOriginalsFasta.log"
    benchmark: "benchmark/D07_GetResolvedOriginalsFasta.tsv"
    shell:"""
        echo "##### D07_GetResolvedOriginalsFasta" > {log}
        {params.workflowDir}/scripts/D07_RenameFastaWithGenes.py {input.fa_unnamed} {input.bed} 1> {output.fa} 2>> {log}
    """
