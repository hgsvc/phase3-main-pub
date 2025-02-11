# Flow of this smk file:
# - E01 Map resolved originals to asm (out: paf). And add Haplotype tag.
# - E02 Needleman Wunch alignment of copy hits to original (calc cigar/matches/mismatches/etc)
# - E02 Calculate resolved copy identity/accuracy
# - E03 Remove hits that intersect with Flagger called regions
# - E04 Filter by identity/accuracy
# - E05 Annotate Original and sort by gene then by hit loc
# - E06 Filter by Original's hit length margin
# - E07 Locate exons
# - E08 Filter by percent alignment of gene model to asm for each copy.
# - E09 Overlapping Nework Filter (to consolidate overlapping genes)
# - E10 paf to bed


rule E01_GetResolvedCopiesPaf:
    input:
        asm="results/A0U_asm.fasta",
        fa="results/D07_resolved_originals.fasta"
    output:
        paf_tmp=temp("results/E01_resolved_copies.tmp.paf"),
        paf="results/E01_resolved_copies.paf"
    resources:
        mem_mb=cluster_mem_mb_xlarge,
        cpus_per_task=cluster_cpus_per_task_large,
        runtime=config["cluster_runtime_long"]
    conda: "../envs/sda2.main.yml"
    log: "logs/E01_GetResolvedCopiesPaf.log"
    benchmark: "benchmark/E01_GetResolvedCopiesPaf.tsv"
    shell:"""
    {{
        echo "##### E01_GetResolvedCopiesPaf" > {log}
        echo "### Minimap originals to asm" >> {log}
        minimap2 -x asm20 -p 0.2 -N 100 -m 10 -E2,0 -s 10 -t {resources.cpus_per_task} {input.asm} {input.fa} 1> {output.paf_tmp}
        echo "### Add haplotype tag"
        cat {output.paf_tmp} | \
            awk 'BEGIN {{OFS="\\t"}} \
                {{if ($6~/haplotype1/) \
                    {{print $0,"hp:Z:haplotype1"}} \
                else if ($6~/haplotype2/) \
                    {{print $0,"hp:Z:haplotype2"}} \
                else \
                    {{print $0}} }}' 1> {output.paf}
    }} 2>> {log}
    """

rule E02_GetResolvedCopyIdentities:
    input:
        paf="results/E01_resolved_copies.paf",
        asm="results/A0U_asm.fasta"
    output:
        pafc="results/E02_mapped_resolved_originals.pafxc",
        pafx="results/E02_mapped_resolved_originals.pafx"
    params:
        workflowDir=workflow.basedir
    resources:
        mem_mb=cluster_mem_mb_medium,
        cpus_per_task=cluster_cpus_per_task_large,
        runtime=config["cluster_runtime_short"],
        tmpdir=tmpDir
    retries: 2
    conda: "../envs/sda2.main.yml"
    log: "logs/E02_GetResolvedCopyIdentities.log"
    benchmark: "benchmark/E02_GetResolvedCopyIdentities.tsv"
    shell:"""
    {{
        echo "##### E02_GetResolvedCopyIdentities" > {log}
        echo "### Align and calculate hit identities" >> {log}
        export OPENBLAS_NUM_THREADS=1        
        cat {input.paf} | xargs -P {resources.cpus_per_task} -I % bash -c ' \
            tmp_dir=`mktemp -d -p {resources.tmpdir} tmp.getIdent.$$.XXXXXX`; \
            echo "$@" | tr "\\t" " " > "$tmp_dir"/line.paf; \
            {params.workflowDir}/scripts/E02_CalcPafIdentity.py {input.asm} "$tmp_dir" "$tmp_dir"/line.paf 1>> {output.pafc}; \
            rm -rf "$tmp_dir"; ' _ %
        cat {output.pafc} | cut -f 1-21 1> {output.pafx} # removes cigar string
    }} 2>> {log}
    """

rule E03_FilterPafByBed:
    input:
        paf="results/E02_mapped_resolved_originals.pafx"
    output:
        bed_draft=temp("results/E03_to_filter_on.bed"),
        filt="results/E03_mapped_resolved_originals_filt.pafx"
    params:
        workflowDir=workflow.basedir,
        bed=config["bed_for_filtering_results"],
        flagger_components=expand("{base}",base=config["flagger_components_to_filter"])
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/E03_FilterPafByBed.log"
    benchmark: "benchmark/E03_FilterPafByBed.tsv"
    shell:"""
    {{
        echo "##### E03_FilterPafByBed" > {log}
        if [[ -f "{params.bed}" && -n "{params.flagger_components}" ]]
        then
            echo "### Run Flagger Filter" >> {log}
            cat {params.bed} | \
                awk 'BEGIN {{OFS="\\t"; data_lines_started="false"}} \
                (data_lines_started=="true" || !/^track/) \
                    {{data_lines_started="true"; \
                    print}}' | \
                awk -v flagger_components="{params.flagger_components}" \
                    'BEGIN \
                        {{OFS="\\t"; \
                        split(flagger_components,fcs_dict," "); \
                        for (i in fcs_dict) fcs[fcs_dict[i]]=""}} \
                    ($4 in fcs) \
                        {{print $1,$2,$3,$4}}' 1> {output.bed_draft}
        elif [[ -f "{params.bed}" && -z "{params.flagger_components}" ]]
        then
            cut -f1-4 {params.bed} 1> {output.bed_draft}
        fi

        if [[ -f "{params.bed}" ]]
        then
            echo "### Run Filter" >> {log}
            cat {input.paf} | \
                awk 'BEGIN {{OFS="\\t"}} {{print $6,$8,$9,$0}}' | \
                bedtools intersect -v -a stdin -b {output.bed_draft} | \
                cut --complement -f1-3 1> {output.filt}
        else
            echo "### Skip filter and generate symbolic link instead." >> {log}
            echo "# Bed to filter by not provided." >> {log}
            ln -s {params.workflowDir}/../{input.paf} {params.workflowDir}/../{output.filt}
            touch {output.bed_draft}
        fi
    }} 2>> {log}
    """

rule E04_FilterLowIdentityPaf:
    input:
        pafx="results/E03_mapped_resolved_originals_filt.pafx"
    output:
        filt="results/E04_mapped_resolved_originals_filtered_by_identity.pafx"
    params:
        min_copy_identity=config["min_copy_identity"]
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/E04_FilterLowIdentityPaf.log"
    benchmark: "benchmark/E04_FilterLowIdentityPaf.tsv"
    shell:"""
    {{
        echo "##### E04_FilterLowIdentityPaf" > {log}
        cat {input.pafx} | \
            awk -v minId={params.min_copy_identity} ' \
                BEGIN \
                    {{OFS="\\t"}} \
                ($20 >= minId) \
                    {{print $0}}' | \
            sort -k1,1 -k6,6 -k8,8n -k9,9n 1> {output.filt}
    }} 2>> {log}
    """

rule E05_AnnotateOriginal:
    input:
        pafx="results/E04_mapped_resolved_originals_filtered_by_identity.pafx"
    output:
        pafxAn="results/E05_mapped_resolved_originals_filtered_ogAnnotated.pafx"
    params:
        marginToCallOriginal=100
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/E05_AnnotateOriginal.log"
    shell:"""
    {{
        echo "##### E05_AnnotateOriginal" > {log}
        cat {input.pafx} | \
            tr '/' '\\t' | \
            awk -v margin={params.marginToCallOriginal} ' \
                function abs(v) {{v += 0; return v < 0 ? -v : v}} \
                BEGIN \
                    {{OFS="\\t"}} \
                {{split($2,og,":"); \
                split(og[2],ogLocs,"-"); \
                hit_chr=$7; hit_start=$9; hit_end=$10; \
                og_chr=og[1]; og_start=ogLocs[1]; og_end=ogLocs[2]; \
                if (og_chr==hit_chr && abs(og_start-hit_start)<margin && abs(og_end-hit_end)<margin) \
                    {{type="Original"}} \
                else \
                    {{type="Copy"}} \
                print $1"/"$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,type}}' 1> {output.pafxAn}
    }} 2>> {log}
    """

# Remove copies that deviate from its original copy's length by 0.1 of the original copy's length.
rule E06_FiltByOriginalLengthMargin:
    input:
        pafx="results/E05_mapped_resolved_originals_filtered_ogAnnotated.pafx"
    output:
        filt="results/E06_mapped_resolved_originals_filtered_ogLength.pafx"
    params:
        maxLengthMargin=config["max_length_margin"]
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/E06_FiltByOriginalLengthMargin.log"
    shell:"""
    {{
        echo "##### E06_FiltByOriginalLengthMargin" > {log}
        cat {input.pafx} | \
            awk -v maxLengthMargin={params.maxLengthMargin} ' \
                BEGIN \
                    {{OFS="\\t"; \
                    highMarg=1+maxLengthMargin; \
                    lowMarg=1-maxLengthMargin}} \
                {{split($1,gene_name,"/"); \
                split(gene_name[2],og,":"); \
                split(og[2],ogLocs,"-"); \
                og_start=ogLocs[1]; og_end=ogLocs[2]; }} \
                (($9-$8)<highMarg*(og_end-og_start) && ($9-$8)>lowMarg*(og_end-og_start)) \
                    {{print $0}}' 1> {output.filt}
    }} 2>> {log}
    """

# Locate Exons in gene copies and filter out gene copies that contain <50% of their gene model
rule E07_LocateExons:
    input:
        pafx="results/E06_mapped_resolved_originals_filtered_ogLength.pafx",
        gm="results/C03_gene_model_filt.fasta",
        asm="results/A0U_asm.fasta"
    output:
        gmidx="results/C03_gene_model_filt.fasta.fai",
        pafxe="results/E07_mapped_resolved_originals_wExons.pafxe"
    params:
        workflowDir=workflow.basedir
    resources:
        mem_mb=cluster_mem_mb_medium,
        cpus_per_task=cluster_cpus_per_task_large,
        runtime=config["cluster_runtime_long"],
        tmpdir=tmpDir
    retries: 2
    conda: "../envs/sda2.main.yml"
    log: "logs/E07_LocateExons.log"
    benchmark: "benchmark/E07_LocateExons.tsv"
    shell:"""
    {{
        echo "##### E07_LocateExons" > {log}
        echo "### Determine Node Variables" >> {log}
        mem_per_cpu="$(echo "{resources.mem_mb}/1.5/{resources.cpus_per_task}" | bc)"
        echo "Memory per cpu: $mem_per_cpu" >> {log}

        echo "### Index Gene Model Fasta" >> {log}
        samtools faidx {input.gm}

        > {resources.tmpdir}/tmp.multiLineTest.err # TODO Delete me

        echo "### Locate Exons" >> {log}
        echo "## Create Locate Function" >> {log}
        locExons () {{ # input = single paf line
            pafLine=`echo "$*" | tr "\\t" " "`

            gene=`echo "$pafLine" | tr "/" "\\t" | awk '(NR==1) {{print $1}}'`
            geneClean=`echo "$pafLine" | tr "/" "\\t" | tr "|" "-" | awk '(NR==1) {{print $1}}'`
            geneLoc=`echo "$pafLine" | tr " " "\\t" | awk '(NR==1) {{print $6":"$8"-"$9}}'`

            tmp_dir=`mktemp -d -p {resources.tmpdir} tmp.getIdent.$geneClean.$$.XXXXXX`

            samtools faidx {input.gm} "$gene" > "$tmp_dir"/gene_model.fasta
            samtools faidx {input.asm} "$geneLoc" > "$tmp_dir"/gene_copy.fasta
            
            minimap2 -x splice -a -t {resources.cpus_per_task} "$tmp_dir"/gene_copy.fasta "$tmp_dir"/gene_model.fasta 1> "$tmp_dir"/gene_copy.sam

            samtools view -b -F 2316 "$tmp_dir"/gene_copy.sam > "$tmp_dir"/gene_copy.bam # -F 4,8,256,2048
            samtools index "$tmp_dir"/gene_copy.bam
            {params.workflowDir}/scripts/D02_AnnotateMappedLength.py "$tmp_dir"/gene_copy.bam | \
                samtools view -b -o "$tmp_dir"/gene_copy_filt.bam

            if [ -f "$tmp_dir"/gene_copy_filt.bam ]; then
                bedtools bamtobed -bed12 -i "$tmp_dir"/gene_copy_filt.bam > "$tmp_dir"/gene_copy.bed
                samtools view "$tmp_dir"/gene_copy_filt.bam | \
                    awk 'BEGIN {{OFS="\\t"}} \
                        {{pa=""; \
                        for (i=12; i<=NF; ++i) \
                            {{if ($i ~ "^pa:f:") \
                                {{pa = substr($i,6,length($i)-5)}} }} \
                        print pa}}' | \
                    paste "$tmp_dir"/gene_copy.bed /dev/stdin/ > "$tmp_dir"/gene_copy.bed13

                head -1 "$tmp_dir"/gene_copy.bed13 | awk -v funcInputPaf="$pafLine" 'BEGIN {{OFS=" "}} {{print funcInputPaf,$11,$12,$13}}' | tr ' ' '\\t' >> {output.pafxe}
                
                cat "$tmp_dir"/gene_copy.bed13 | tail -n+2 >> {resources.tmpdir}/tmp.multiLineTest.err
            fi

            rm -rf "$tmp_dir"
        }}
        export -f locExons

        echo "## Run Locate Function" >> {log}
        cat {input.pafx} | tr "\\t" " " | \
            xargs -P {resources.cpus_per_task} -I % bash -c ' \
                locExons "$@" ' _ %
    }} 2>> {log}
    """

rule E08_FilterByModelMappedLength:
    input:
        pafxe="results/E07_mapped_resolved_originals_wExons.pafxe"
    output:
        filt="results/E08_mapped_resolved_originals_filtered_modelMappedLength.pafxe",
        igv_bed="results/E08_mapped_resolved_originals_wExons.igv.bed"
    params:
        min_gm_alignment=config["min_gene_model_alignment"],
        workflowDir=workflow.basedir
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/E08_FilterByModelMappedLength.log"
    benchmark: "benchmark/E08_FilterByModelMappedLength.tsv"
    shell:"""
    {{
        echo "##### E08_FilterByModelMappedLength" > {log}
        cat {input.pafxe} | \
            awk 'BEGIN {{OFS="\\t"}} \
                ($25>={params.min_gm_alignment}) \
                    {{print}}' > {output.filt}

        echo "### IGV BED File" >> {log}
        cat {output.filt} | \
            awk 'BEGIN {{OFS="\\t"}} \
                {{split($22,exonSizes,","); \
                print $6,$8,$9,$1,int($20*1000),$5,$8,$9,"0,0,255",length(exonSizes),$23,$24}}' 1> {output.igv_bed}
    }} 2>> {log}
    """

rule E09_GroupIsoformsByExonOverlap:
    input:
        pafxe="results/E08_mapped_resolved_originals_filtered_modelMappedLength.pafxe"
    output:
        marked="results/E09_mapped_resolved_originals_annotated_network.pafxe",
        filt="results/E09_mapped_resolved_originals_filtered_network.pafxe",
        coms="results/E09_isoform_communities_groupedByExonOverlap.tsv"
    params:
        uncharacterized_prefix=config["uncharacterized_gene_name_prefix"],
        workflowDir=workflow.basedir
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/E09_GroupIsoformsByExonOverlap.log"
    benchmark: "benchmark/E09_GroupIsoformsByExonOverlap.tsv"
    shell:"""
    {{
        echo "##### E09_GroupIsoformsByExonOverlap" > {log}
        echo "### Identify Prefix of Uncharacterized Genes" >> {log}
        prefix="{params.uncharacterized_prefix}"
        prefix_len=${{#prefix}}
        if [[ $prefix_len -gt 0 ]]
        then
            echo "# Deprioritizing of gene names with prefix \'{params.uncharacterized_prefix}\' enabled when picking a representative gene name within a community as recognized bythe exon network filter." >> {log}
            unchar_filt_flag="-u {params.uncharacterized_prefix}"
        else
            echo "# Deprioritization of uncharacterized genes when picking a representative gene name is disabled." >> {log}
            unchar_filt_flag=""
        fi

        echo "### Annotate non-representative overlapping genes" >> {log}
        {params.workflowDir}/scripts/E09_NetworkFilter.py {input.pafxe} {output.coms} "$unchar_filt_flag" 1> {output.marked}
        cat {output.marked} | \
            awk 'BEGIN {{OFS="\\t"}} ($26=="yes") {{print $0}}' | \
            cut -f1-25 1> {output.filt}
    }} 2>> {log}
    """

rule E10_FinalResolvedCopiesBed:
    input:
        pafxe="results/E08_mapped_resolved_originals_filtered_modelMappedLength.pafxe",
        pafx="results/E09_mapped_resolved_originals_filtered_network.pafxe"
    output:
        bed="results/E10_resolved_copies.bed", # unfiltered by E09 Network Filter Results
        bed12="results/E10_resolved_copies.igv.bed" # filtered according to E09 Network Filter Results
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/E10_FinalResolvedCopiesBed.log"
    shell:"""
    {{
        echo "##### E10_FinalResolvedCopiesBed" > {log}

        cat {input.pafxe} | \
            awk 'BEGIN {{OFS="\\t"}} \
                {{if ($5=="+") \
                    {{strand=0}} \
                else \
                    {{strand=1}} \
                split($1,gene_name,"/"); \
                split(gene_name[2],og,":"); \
                split(og[2],ogLocs,"-"); \
                og_chr=og[1]; og_start=ogLocs[1]; og_end=ogLocs[2]; \
            print $6,$8,$9,gene_name[1],og_chr,og_start,og_end,strand,$19,$20,$21,$22,$23,$24,$25}}' | \
            sort -k1,1 -k2,2n -k3,3n -k4,4 1> {output.bed}
        # Out format: chr,start,end,gene,original_chr,original_start,original_end,strand,p_identity,p_accuracy,exons_sizes,exon_starts,gm_alignment
        # sorted by chrom, start, end, then gene

        echo "### IGV BED File" >> {log}
        cat {input.pafx} | \
            awk 'BEGIN {{OFS="\\t"}} \
                {{split($23,exonSizes,","); \
                print $6,$8,$9,$1,int($20*1000),$5,$8,$9,"0,255,0",length(exonSizes),$23,$24}}' 1> {output.bed12}
    }} 2>> {log}
    """
