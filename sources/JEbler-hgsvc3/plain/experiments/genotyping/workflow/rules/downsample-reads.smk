downsampling_reads = {}
for line in open(READS, 'r'):
        if line.startswith('#'):
                continue
        fields = line.strip().split()
        downsampling_reads[fields[1]] = fields[7] 


# align the reads to the reference to estimate the coverage of the full data set
# need to do it per callset because the reference might be different
rule align_reads:
	input:
		reads = lambda wildcards: downsampling_reads[wildcards.sample],
		fasta = lambda wildcards: CALLSETS[wildcards.callset]['reference'] 
	output:
		bam = "{results}/downsampling/{callset}/{coverage}/aligned/{sample}_full.bam",
		bai = "{results}/downsampling/{callset}/{coverage}/aligned/{sample}_full.bam.bai"
	conda:
		"../envs/downsampling.yml"
	resources:
		mem_total_mb = 60000,
		runtime_hrs = 25,
		runtime_min = 1
	log:
		"{results}/downsampling/{callset}/{coverage}/aligned/{sample}_full.log"
	shell:
		'(/usr/bin/time -v bwa mem -t {threads} -M {input.fasta} -R "@RG\\tID:{wildcards.sample}\\tLB:lib1\\tPL:illumina\\tPU:unit1\\tSM:{wildcards.sample}" {input.reads} | samtools view -bS | samtools sort -o {output} - ) &> {log}'


# estimate the coverage of the aligned data
rule compute_bam_coverage:
	input:
		full_cov_bam = "{results}/downsampling/{callset}/{coverage}/aligned/{sample}_full.bam",
		full_cov_bai = "{results}/downsampling/{callset}/{coverage}/aligned/{sample}_full.bam.bai"
	output:
		"{results}/downsampling/{callset}/{coverage}/aligned/{sample}_full.cov"
	conda:
		'../env/downsampling.yml'
	resources:
		mem_total_mb = 10000,
		runtime_hrs = 5,
		runtime_min = 1
	log:
		"{results}/downsampling/{callset}/{coverage}/aligned/{sample}_full_cov.log"
	shell:
		"bash workflow/scripts/compute-coverage.sh {input.full_cov_bam} {output} &> {log}"


# downsample fastq to desired coverage
rule downsample_reads:
	input:
		lambda wildcards: downsampling_reads[wildcards.sample]
	output:
		"{results}/downsampling/{callset}/{coverage}/{sample}_{coverage}.fa.gz"
	conda:
		"../envs/downsampling.yml"
	resources:
		mem_total_mb = 20000,
		runtime_hrs = 5,
		runtime_min = 1
	log:
		"{results}/downsampling/{callset}/{coverage}/{sample}_{coverage}.log"
	shell:
		"bash workflow/scripts/downsample-fasta.sh {input.coverage} {wildcards.fraction} {input} {output} &> {log}"
