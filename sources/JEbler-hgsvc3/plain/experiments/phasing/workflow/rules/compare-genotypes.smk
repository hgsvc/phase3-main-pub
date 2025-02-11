


rule extract_sample_phasing_with_absent:
	"""
	Extract phasing of a single sample. This is done to reduce the
	time/memory of whatshap compare.
	"""
	input:
		"{results}/{version}/phased_{version}.vcf.gz"
	output:
		vcf = "{results}/concordance/{version}/phased_{version}_{sample}_all.vcf.gz",
		tbi = "{results}/concordance/{version}/phased_{version}_{sample}_all.vcf.gz.tbi"
	conda:
		"../envs/shapeit.yaml"
	wildcard_constraints:
		version = "shapeit"
	resources:
		mem_total_mb = 50000,
		runtime_hrs = 2
	params:
		chrom = ','.join([c for c in MAPS.keys() if (not 'X' in c) and (not 'Y' in c)])
	threads:
		10
	shell:
		"""
		bcftools view --threads {threads} --samples {wildcards.sample} --regions {params.chrom} {input} | bgzip > {output.vcf}
		tabix -p vcf {output.vcf}
		"""


rule extract_sample_truthset_with_absent:
	"""
	Extract ground truth phasing of a single sample. Done to reduce
	the time/memory of whatshap compare.
	"""
	input:
		lambda wildcards: TRUTHSETS[wildcards.truthset]["vcf"]
	output:
		vcf = "{results}/concordance/{truthset}_{sample}_all.vcf.gz",
		tbi = "{results}/concordance/{truthset}_{sample}_all.vcf.gz.tbi"
	conda:
		"../envs/shapeit.yaml"
	wildcard_constraints:
		truthsets = "|".join(TRUTHSETS.keys()),
		sample = "|".join(EVAL_SAMPLES)
	resources:
		mem_total_mb = 50000,
		runtime_hrs = 2
	threads:
		10
	params:
		chrom = ','.join([c for c in MAPS.keys() if (not 'X' in c) and (not 'Y' in c)])
	shell:
		"""
		bcftools view --samples {wildcards.sample} --threads {threads} --regions {params.chrom} {input} | bgzip > {output.vcf}
		tabix -p vcf {output.vcf}
		"""
	


rule compare_genotypes:
	input:
		truthset = "{results}/concordance/{truthset}_{sample}_all.vcf.gz",
		phasing = "{results}/concordance/{version}/phased_{version}_{sample}_all.vcf.gz"
	output:
		"{results}/concordance/{truthset}_{version}_{sample}.txt"
	resources:
		mem_total_mb = 40000,
		runtime_hrs = 0,
		runtime_min = 40
	conda:
		"../envs/shapeit.yaml"
	shell:
		"""
		python3 workflow/scripts/genotype-evaluation-intersection.py {input.truthset} {input.phasing} > {output}
		"""


rule plot_genotype_comparison:
	input:
		lambda wildcards: expand("{{results}}/concordance/{{truthset}}_{{version}}_{sample}.txt", sample = TRUTHSETS[wildcards.truthset]["evaluation_samples"])
	output:
		pdf = "{results}/concordance/{truthset}_{version}.pdf",
		tsv = "{results}/concordance/{truthset}_{version}.tsv"
	conda:
		"../envs/plotting.yaml"
	shell:
		"""
		ls {input} | python3 workflow/scripts/plot-concordances.py {output.pdf} {wildcards.truthset} > {output.tsv}
		"""
