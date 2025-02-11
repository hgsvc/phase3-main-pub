rule B00A_GenerateFilesWithoutUsingBAMs:
    input:
        fai="results/A0U_{hap}_asm.fasta.fai"
    output:
        cov="results/B01_{hap}_hmm/B01_cov_bins.bed.gz",
        vcfDups="results/B02_{hap}_copy_number.bed.gz"
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/B00A_GenerateFilesWithoutUsingBAMs.{hap}.log"
    benchmark: "benchmark/B00A_GenerateFilesWithoutUsingBAMs.{hap}.tsv"
    shell:"""
    {{
        echo "##### B00A_GenerateFilesWithoutUsingBAMs - {wildcards.hap}" > {log}
        echo "### Generate B01_cov_bins.bed.gz" >> {log}
        cat {input.fai} | \
            sort -k1,1 | \
            awk 'BEGIN {{OFS="\\t"}} \
                {{for (i=0; i+100 <= $2; i+=100) \
                    {{print $1,i,i+100,100;}} }}' | \
            bgzip -c 1> {output.cov}
            # Out format: chr,start,end,read_depth
        
        echo "### Generate B02_copy_number.bed.gz" >> {log}
        cat {input.fai} | \
            awk 'BEGIN {{OFS="\\t"}} \
                {{print $1,1,$2,1,100}}' | \
            bgzip -c 1> {output.vcfDups}
            # Out format: chr,start,end,copy_num,read_depth
    }} 2>> {log}
    """

rule B00B_GenerateFilesWithoutUsingBAMs:
    output:
        warning="WARNING_NO_REAL_BAMS_USED.txt",
        mean="results/B03_asm_mean_cov.txt",
        mean2="results/B06_coding_asm_mean_cov.txt",
        hapChrs="results/B05_haploid_chrs.txt"
    params:
        autodetect_haploid_chrs_flag=config["flag_autodetect_haploid_chrs"],
        hapChrs=expand("{base}",base=config["haploid_chrs"])
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/B00B_GenerateFilesWithoutUsingBAMs.log"
    benchmark: "benchmark/B00B_GenerateFilesWithoutUsingBAMs.tsv"
    shell:"""
    {{
        echo "##### B00B_GenerateFilesWithoutUsingBAMs" > {log}
        echo "### Generate Warning Notices" >> {log}
        echo "*** WARNING ***"
        echo "No read bam paths provided!!"
        echo "This pipeline will run assuming depth across the assembly is 100 bp."
        echo "*** WARNING ***" >> {log}
        echo "No read bam paths provided!!" >> {log}
        echo "This pipeline will run assuming depth across the assembly is 100 bp." >> {log}
        echo "*** WARNING ***" >> {output.warning}
        echo "No read bam paths provided!!" >> {output.warning}
        echo "This pipeline will run assuming depth across the assembly is 100 bp." >> {output.warning}

        echo "### Generate B03_asm_mean_cov.txt" >> {log}
        echo "100" > {output.mean}
        echo "100" > {output.mean2}

        echo "### Generate B05_haploid_chrs.txt" >> {log}
        if [ {params.autodetect_haploid_chrs_flag} = "True" ]
        then
            echo "# Haploid chromosome autodetection set." >> {log}
            echo "Given no read data, all chromosomes are assumed to be diploid." >> {log}
            echo "Given chromosome autodetection set without any read data, all chromosomes are assumed to be diploid."
            echo "Given chromosome autodetection set without any read data, all chromosomes are assumed to be diploid." >> {output.warning}
            > {output.hapChrs}
        else
            echo "# Using supplied haploid chromosomes" >> {log}
            echo {params.hapChrs} | tr ' ' '\n' 1> {output.hapChrs}
        fi
    }} 2>> {log}
    """
