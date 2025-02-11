rule concat_vcfs:
	"""
	Combine the phased per-chromosome VCFs into a single one.
	"""
	input:
		expand("{{results}}/{{version}}/phased_{{version}}_{chrom}.bcf", chrom = [c for c in MAPS.keys()])
	output:
		"{results}/{version}/phased_{version}.vcf.gz"
	conda:
		"../envs/shapeit.yaml"
	wildcard_constraints:
		version = "shapeit|shapeit-scaffold"
	benchmark:
		"{results}/{version}/phased_{version}_benchmark.txt"
	threads: 24
        resources:
		mem_total_mb = 100000,
		runtime_hrs = 5,
		runtime_min = 1
	log:
		"{results}/{version}/phased_{version}.log"
	shell:
		"""
		bcftools concat -o {output} -O z --threads {threads} {input} &> {log}
		tabix -p vcf {output}
		"""


rule extract_sample_phasing:
	"""
	Extract phasing of a single sample. This is done to reduce the
	time/memory of whatshap compare.
	"""
	input:
		"{results}/{version}/phased_{version}.vcf.gz"
	output:
		vcf = temp("{results}/evaluation/{version}/phased_{version}_{sample}_all.vcf.gz"),
		tbi = temp("{results}/evaluation/{version}/phased_{version}_{sample}_all.vcf.gz.tbi")
	conda:
		"../envs/shapeit.yaml"
	resources:
		mem_total_mb = 50000,
		runtime_hrs = 2
	params:
		chrom = ','.join([c for c in MAPS.keys() if (not 'X' in c) and (not 'Y' in c)])
	threads:
		10
	shell:
		"""
		bcftools view --threads {threads} --samples {wildcards.sample} --regions {params.chrom} {input} | bcftools view --min-ac 1 | bgzip > {output.vcf}
		tabix -p vcf {output.vcf}
		"""


rule extract_sample_truthset:
	"""
	Extract ground truth phasing of a single sample. Done to reduce
	the time/memory of whatshap compare.
	"""
	input:
		lambda wildcards: TRUTHSETS[wildcards.truthset]["vcf"]
	output:
		vcf = temp("{results}/evaluation/{truthset}_{sample}_all.vcf.gz"),
		tbi = temp("{results}/evaluation/{truthset}_{sample}_all.vcf.gz.tbi")
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
		bcftools view --samples {wildcards.sample} --threads {threads} --regions {params.chrom} {input} | bcftools view --min-ac 1 | bgzip > {output.vcf}
		tabix -p vcf {output.vcf}
		"""
		

rule prepare_ps_header:
	output:
		"{results}/evaluation/header_ps.txt"
	shell:
		"""
		echo "##INFO=<ID=PS,Number=1,Type=Integer,Description=\\"Phase set.\\">" > {output}
		"""



rule extract_region:
	"""
	Extract a region and add different PS tag for each BED interval.
	In this way, switch errors between regions will not be penalized.
	"""
	input:
		vcf = "{results}/evaluation/{filename}_all.vcf.gz",
		bed = lambda wildcards: PHASING_REGIONS[wildcards.region],
		header = "{results}/evaluation/header_ps.txt"
	output:
		vcf = "{results}/evaluation/{filename}_{region}.vcf.gz",
		tbi = "{results}/evaluation/{filename}_{region}.vcf.gz.tbi",
		bed = temp("{results}/evaluation/{filename}_{region}.bed.gz")
	conda:
		"../envs/shapeit.yaml"
	wildcard_constraints:
		region = "|".join(PHASING_REGIONS.keys())
	shell:
		"""
		bedtools merge -i {input.bed} | python3 workflow/scripts/prepare-bed.py | bgzip > {output.bed}
		tabix -p bed {output.bed}
		bcftools view -R {output.bed} {input.vcf} |  bcftools annotate -a {output.bed} -c CHROM,FROM,TO,PS -h {input.header} | python3 workflow/scripts/transfer-ps.py | bgzip > {output.vcf}
		tabix -p vcf {output.vcf}
		"""
	


rule evaluate:
	"""
	Evaluate Phasing based on truth sets.
	"""
	input:
		truthsets = "{results}/evaluation/{truthset}_{sample}_{region}.vcf.gz",
		truthsets_tbi = "{results}/evaluation/{truthset}_{sample}_{region}.vcf.gz.tbi",
		phasing = "{results}/evaluation/{version}/phased_{version}_{sample}_{region}.vcf.gz",
		phasing_tbi = "{results}/evaluation/{version}/phased_{version}_{sample}_{region}.vcf.gz.tbi"
	output:
		pair = "{results}/evaluation/{version}/evaluation_{version}_{sample}_{truthset}_{vartype}_{region}_pair.tsv",
		multi = "{results}/evaluation/{version}/evaluation_{version}_{sample}_{truthset}_{vartype}_{region}_multi.tsv"
	conda:
		"../envs/whatshap.yaml"
	wildcard_constraints:
		truthset = "|".join(TRUTHSETS.keys()),
		sample = "|".join(EVAL_SAMPLES),
		vartype = "snps|all"
	resources:
		mem_total_mb = 100000,
		runtime_hrs = 3
	params:
		names = lambda wildcards: wildcards.version + "," + wildcards.truthset,
		variants = lambda wildcards: "--only-snvs" if wildcards.vartype == "snps" else ""
	log:
		"{results}/evaluation/{version}/evaluation_{version}_{sample}_{truthset}_{vartype}_{region}.log"
	shell:
		"""
		whatshap compare {input.phasing} {input.truthsets} --sample {wildcards.sample} --tsv-pairwise {output.pair} --tsv-multiway {output.multi} --names {params.names} {params.variants} &> {log}
		"""


rule plot_phasing_results:
	"""
	Plot the switch error rates.
	"""
	input:
		tsv_files = lambda wildcards: expand("{{results}}/evaluation/{{version}}/evaluation_{{version}}_{sample}_{{truthset}}_{{vartype}}_{region}_pair.tsv", sample = TRUTHSETS[wildcards.truthset]["evaluation_samples"], region=["all"] + [r for r in PHASING_REGIONS.keys()]),
		trios = FAM 
	output:
		"{results}/evaluation/{version}/evaluation_{version}_{truthset}_{vartype}.pdf"
	log:
		"{results}/evaluation/{version}/evaluation_{version}_{truthset}_{vartype}.log"
	wildcard_constraints:
		truthset='|'.join([c for c in TRUTHSETS.keys()]),
		vartype = "snps|all"
	conda:
		"../envs/plotting.yaml"
	shell:
		"""
		python3 workflow/scripts/plot-phasing-results.py -tsvfiles {input.tsv_files} -truthsetname {wildcards.truthset} -trios {input.trios} -outname {output} &> {log}
		"""


rule multiway_evaluation:
	input:
		truthsets = expand("{{results}}/evaluation/{truthset}_{{sample}}_{{region}}.vcf.gz", truthset = [k for k in config["truthsets"].keys()]),
		truthsets_tbi = expand("{{results}}/evaluation/{truthset}_{{sample}}_{{region}}.vcf.gz.tbi", truthset = [k for k in config["truthsets"].keys()]),
		phasing = "{results}/evaluation/{version}/phased_{version}_{sample}_{region}.vcf.gz",
		phasing_tbi = "{results}/evaluation/{version}/phased_{version}_{sample}_{region}.vcf.gz.tbi"
	output:
		pair = "{results}/evaluation/{version}/evaluation_{version}_{sample}_{vartype}_{region}_multiway_pair.tsv",
		multi = "{results}/evaluation/{version}/evaluation_{version}_{sample}_{vartype}_{region}_multiway_multi.tsv"
	conda:
		"../envs/whatshap.yaml"
	wildcard_constraints:
		sample = "|".join(EVAL_SAMPLES)
	resources:
		mem_total_mb = 50000,
		runtime_hrs = 3
	params:
		names = "shapeit-phasing," + ",".join([k for k in TRUTHSETS.keys()]),
		variants = lambda wildcards: "--only-snvs" if wildcards.vartype == "snps" else ""
	log:
		"{results}/evaluation/evaluation_{version}_{sample}_{vartype}_{region}_multiway.log"
	shell:
		"""
		whatshap compare {input.phasing} {input.truthsets} --sample {wildcards.sample} --tsv-pairwise {output.pair} --tsv-multiway {output.multi} --names {params.names} {params.variants} &> {log}
		"""


rule plot_multiway:
	input:
		tsv=expand("{{results}}/evaluation/{{version}}/evaluation_{{version}}_{sample}_{{vartype}}_{{region}}_multiway_multi.tsv", sample=JOINT_EVAL_SAMPLES),
		trios=FAM
	output:
		"{results}/evaluation/{version}/evaluation_{version}_{vartype}_{region}_multiway.pdf"
	conda:
		"../envs/plotting.yaml"
	shell:
		"""
		python3 workflow/scripts/plot-multiway.py -tsv {input.tsv} -trios {input.trios} -outname {output} -title {wildcards.region}
		"""
