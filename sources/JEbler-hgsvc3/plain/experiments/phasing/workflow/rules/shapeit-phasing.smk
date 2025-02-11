rule shapeit_extract_males:
	input:
		SEX
	output:
		"{results}/haploid-samples.txt"
	shell:
		"awk '$2==1' {input} > {output}"


rule shapeit_extract_chromosome:
	"""
	Extract specific chromosome.
	"""
	input:
		GENOTYPES
	output:
		"{results}/vcf/genotypes_{chrom}.vcf.gz"
	conda:
		"../envs/shapeit.yaml"
	log:
		"{results}/vcf/genotypes_{chrom}.log"
	resources:
		mem_total_mb = 70000,
		runtime_hrs = 5
	threads: 10
	shell:
		"""
		bcftools view -r {wildcards.chrom} {input} --threads {threads} -Oz -o {output} &> {log}
		tabix -p vcf {output}
		"""
	

rule shapeit_set_low_qual_to_missing:
	"""
	Set low quality genotypes to missing.
	"""
	input:
		 lambda wildcards: MERGED_GENOTYPES[wildcards.chrom] if MERGED_GENOTYPES else "{results}/vcf/genotypes_{chrom}.vcf.gz" 
	output:
		temp("{results}/vcf/low-qual-missing_{chrom}.vcf.gz")
	log:
		"{results}/vcf/low-qual-missing_{chrom}.log"
	conda:
		"../envs/shapeit.yaml"
	resources:
		mem_total_mb = 50000,
		runtime_hrs = 5
	shell:
		"""
		bcftools +setGT {input}  -- -t q -n ./. -i 'FMT/GQ<10' 2> {log} | bgzip > {output}
		tabix -p vcf {output}
		"""

def phasing_input(wildcards):
	if GQ_FILTERING:
		return "{results}/vcf/low-qual-missing_{chrom}.vcf.gz".format(results=wildcards.results, chrom=wildcards.chrom)
	elif MERGED_GENOTYPES:
		return MERGED_GENOTYPES[wildcards.chrom]
	else:
		return "{results}/vcf/genotypes_{chrom}.vcf.gz".format(results=wildcards.results, chrom=wildcards.chrom)


rule shapeit_phase_common:
	"""
	Phase a chromosome using shapeit.
	"""
	input:
		vcf = phasing_input,
		fam = FAM,
		map = lambda wildcards: MAPS[wildcards.chrom],
		haploids = "{results}/haploid-samples.txt"
	output:
		"{results}/shapeit/phased_shapeit_{chrom}.bcf"
	log:
		"{results}/shapeit/phased_shapeit_{chrom}.log"
	benchmark:
		"{results}/shapeit/phased_shapeit_{chrom}-benchmark.txt"
	conda:
		"../envs/shapeit.yaml"
	threads: 32
        resources:
		mem_total_mb = 200000, # 100000
		runtime_hrs = 30
	params:
		haploids = lambda wildcards: "--haploids " + "{results}/haploid-samples.txt".format(results = wildcards.results)  if ("X" in wildcards.chrom) or ("Y" in wildcards.chrom) else ""
	shell:
		"""
		SHAPEIT5_phase_common --input {input.vcf} --pedigree {input.fam} --region {wildcards.chrom} {params.haploids} --map {input.map} --output {output} --thread {threads} &> {log}
		"""

