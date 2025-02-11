# Flow of this smk file:
# - F01 Calculate weighted depth over gene w/o hmm
# - F02 Label Gene Copy Nums
# - F03 Calculate depth over gene given hmm
# - F04 Filter by min_depth (default: 0.05)
# - F05 Filter by copy number or collapse presence:
#          - autosomal > 2 copies
#          - sex chr > 1 copy
#          - autosomes + sex chr > 1 copy
#          - resolved duplication > 1 copy (hap_unaware)
#          - resolved duplication > 2 copies (hap_aware)
#          - resolved duplication > 1 copy on either haplotype (hap_aware)
# - F06 Group Isoforms that have any overlap
# - F07 Group Isoforms by label created in E09

rule F01_GetGeneCoverage: # naive depths
    input:
        depths=expand("results/B01_{hap}_hmm/B01_cov_bins.bed.gz", hap=OUTPUT_HAPLOTYPES),
        genes="results/E10_resolved_copies.bed" # presorted
    output:
        bed="results/F01_resolved_copies_cn.bed"
    resources:
        mem_mb=cluster_mem_mb_large,
        cpus_per_task=cluster_cpus_per_task_medium,
        runtime=config["cluster_runtime_long"]
    conda: "../envs/sda2.main.yml"
    log: "logs/F01_GetGeneCoverage.log"
    benchmark: "benchmark/F01_GetGeneCoverage.tsv"
    shell:"""
    {{
        echo "##### F01_GetGeneCoverage" > {log}
        zcat {input.depths} | \
            sort -k1,1 -k2,2n | \
            bedtools intersect -loj -sorted -a {input.genes} -b stdin | \
            bedtools groupby -g 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 -c 19,19 -o mean,stdev 1> {output.bed}
        # Out format: chr,start,end,gene,original_chr,original_start,original_end,strand,haplotype,p_identity,p_accuracy,copy,exons_sizes,exon_starts,gm_alignment,depth_by_traditional(non-vcf),depth_by_traditional(non-vcf)-stdDev
    }} 2>> {log}
    """

# Label Copy Nums and convert depth from bp to fraction of mean asm
rule F02_LabelCopyNums:
    input:
        bed="results/F01_resolved_copies_cn.bed",
        mean="results/B06_coding_asm_mean_cov.txt"
    output:
        bed="results/F02_resolved_copies_cn_CNsLabeled.bed"
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/F02_LabelCopyNums.log"
    benchmark: "benchmark/F02_LabelCopyNums.tsv"
    shell:"""
    {{
        echo "##### F02_LabelCopyNums" > {log}
        cat {input.bed} | \
            awk -v meanCov=$(cat {input.mean}) ' \
                function round(x) {{i=int(x+0.5); \
                                    if (i<1) \
                                        {{return 1}} \
                                    else \
                                        {{return i}}}} \
                BEGIN \
                    {{OFS="\\t"}} \
                {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16/meanCov,$17/meanCov,round($16/meanCov)}}' 1> {output.bed}
    }} 2>> {log}
    """

# TODO verify that "." is 0 and not 1, I think it's 1 but maybe not all the time?
rule F03_GetGeneHmmCoverage: # hmm vcf depths
    input:
        hmm=expand("results/B02_{hap}_copy_number.bed.gz", hap=OUTPUT_HAPLOTYPES),
        genes="results/F02_resolved_copies_cn_CNsLabeled.bed" # presorted
    output:
        hmm_noZero=temp("results/F03_copy_number_filtered.bed"),
        bed="results/F03_resolved_copies_cn2.bed"
    resources:
        mem_mb=cluster_mem_mb_large,
        cpus_per_task=cluster_cpus_per_task_medium,
        runtime=config["cluster_runtime_long"]
    conda: "../envs/sda2.main.yml"
    log: "logs/F03_GetGeneHmmCoverage.log"
    benchmark: "benchmark/F03_GetGeneHmmCoverage.tsv"
    shell:"""
    {{
        echo "##### F03_GetGeneHmmCoverage" > {log}
        zcat {input.hmm} | \
            sort -k1,1 -k2,2n | \
            awk 'BEGIN {{OFS="\\t"}} ($4!=0) {{print}}' 1> {output.hmm_noZero}
        bedtools intersect -loj -a {input.genes} -b {output.hmm_noZero} -f 1 | \
            awk 'BEGIN {{OFS="\\t"}} \
                {{if ($21==".") \
                    {{$21=0}} \
                print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$22}}' 1> {output.bed}
    }} 2>> {log}
    """

rule F04_FilterOutMinDepth:
    input:
        bed="results/F03_resolved_copies_cn2.bed"
    output:
        filt="results/F04_resolved_copies_cn2_minDepthFilt.bed"
    params:
        minDepth=config["min_depth"]
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/F04_FilterOutMinDepth.log"
    benchmark: "benchmark/F04_FilterOutMinDepth.tsv"
    shell:"""
    {{
        echo "##### F04_FilterOutMinDepth" > {log}
            cat {input.bed} | \
            awk -v minDepth={params.minDepth} ' \
                BEGIN \
                    {{OFS="\\t"}} \
                ($16>minDepth) \
                    {{print}}' 1> {output.filt}
    }} 2>> {log}
    """

if hap_aware_mode:
    # - F05H Filter by copy number or collapse presence:
    #          - >= 3 copies (collapsed plus resolved)
    #          - >= 1 collapsed copy
    #          - resolved duplication > 1 copy on either haplotype
    rule F05H_FindDups:
        input:
            bed="results/F04_resolved_copies_cn2_minDepthFilt.bed"
        output:
            dups="results/F05_dups_allFams.bed"
        params:
            workflowDir=workflow.basedir
        localrule: True
        conda: "../envs/sda2.main.yml"
        log: "logs/F05H_FindDups.log"
        benchmark: "benchmark/F05H_FindDups.tsv"
        shell:"""
        {{
            echo "##### F05H_FindDups" > {log}
            echo "Haplotype aware mode activated." >> {log}
            # sort so OG before copies and all genes adjacent
            # traverse and isolate by gene groups
                # if copyNum of any copy > 1 keep gene group
                # if resolved copies present keep gene group

            cat {input.bed} | \
                sort -k4,4 -k12,12r -k5,5 -k6,6n -k7,7n -k1,1 -k2,2n -k3,3n | \
                {params.workflowDir}/scripts/F05_SelectDuplications.py /dev/stdin --phased_input 1> {output.dups}
        }} 2>> {log}
        """
else:
    # - F05 Filter by copy number or collapse presence:
    #          - sex chr > 1 copy
    #          - autosomes + sex chr > 1 copy
    #          - autosomes >= 2 copies
    #          - resolved duplication > 1 copy
    #          - collapsed copy >= 1
    rule F05_FindDups:
        input:
            bed="results/F04_resolved_copies_cn2_minDepthFilt.bed",
            hapChrs="results/B05_haploid_chrs.txt"
        output:
            dups="results/F05_dups_allFams.bed"
        params:
            workflowDir=workflow.basedir
        localrule: True
        conda: "../envs/sda2.main.yml"
        log: "logs/F05_FindDups.log"
        benchmark: "benchmark/F05_FindDups.tsv"
        shell:"""
        {{
            echo "##### F05_FindDups" > {log}
            echo "Haplotype aware mode disabled." >> {log}
            # sort so OG before copies and all genes adjacent
            # traverse and isolate by gene groups
                # if copyNum of any copy > 1 or 2 (depending on if on given haploid chromosome) keep gene group
                # if resolved copies present keep gene group

            cat {input.bed} | \
                sort -k4,4 -k12,12r -k5,5 -k6,6n -k7,7n -k1,1 -k2,2n -k3,3n | \
                {params.workflowDir}/scripts/F05_SelectDuplications.py /dev/stdin --sex_chrs_list_filepath {input.hapChrs} 1> {output.dups}
        }} 2>> {log}
        """

rule F06_GroupIsoformsByAnyOverlap:
    input:
        bed="results/F05_dups_allFams.bed"
    output:
        tsv=temp("results/F06_isoforms_groupedByAnyOverlap_initial.tsv"),
        coms="results/F06_isoform_communities_groupedByAnyOverlap.tsv",
        bed_filt="results/F06_dups_groupedByAnyOverlap.bed"
    params:
        workflowDir=workflow.basedir
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/F06_GroupIsoformsByAnyOverlap.log"
    shell:"""
    {{
        echo "##### F06_GroupIsoformsByAnyOverlap" > {log}
        cat {input.bed} | \
            sort -k1,1 -k2,2n -k3,3n | \
            awk 'BEGIN {{OFS="\\t"; chrm=""; groupEnd=0; isos=""; isos_c=""}} \
                (NR==1) \
                    {{chrm=$1; groupEnd=$3}} \
                ($1!=chrm || $2>groupEnd) \
                    {{print isos,isos_c; \
                    chrm=$1; groupEnd=$3; isos=""; isos_c=""}} \
                {{isos=isos $4","; isos_c=isos_c $4"/"$1":"$2"-"$3","; \
                if ($3>groupEnd) \
                    {{groupEnd=$3}}}} \
                END \
                    {{print isos,isos_c}}' | \
            sed 's/,\\t/\\t/g' | \
            sed 's/,$//g' > {output.tsv}

        {params.workflowDir}/scripts/F06_mergeIsoFamsAndPickConsensus.py {output.tsv} {output.coms} {input.bed} > {output.bed_filt}
    }} 2>> {log}
    """

rule F07_PickRepresentativeForIsoformsGroupedByExonOverlap:
    input:
        bed="results/F05_dups_allFams.bed",
        iso="results/E09_isoform_communities_groupedByExonOverlap.tsv"
    output:
        iso_simp=temp("results/F07_isoform_communities_groupedByExonOverlap_initial_simplified.tsv"),
        coms="results/F07_isoform_communities_groupedByExonOverlap.tsv",
        reps="results/F07_dups_groupedByExonOverlap.bed"
    params:
        workflowDir=workflow.basedir
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/F07_PickRepresentativeForIsoformsGroupedByExonOverlap.log"
    shell:"""
    {{
        echo "##### F07_PickRepresentativeForIsoformsGroupedByExonOverlap" > {log}
        cat {input.iso} | \
            tail -n+2 | \
            cut -f2 | \
            awk 'BEGIN {{OFS="\\t"}} \
                {{n=patsplit($1,ogs,"/[^,]*[,]?",isos); \
                j=isos[0]; \
                for (i=1; i<n; i++) \
                    {{j=j "," isos[i]}}; \
                print j,$1}}' 1> {output.iso_simp}

        {params.workflowDir}/scripts/F06_mergeIsoFamsAndPickConsensus.py {output.iso_simp} {output.coms} {input.bed} 1> {output.reps}
    }} 2>> {log}
    """

