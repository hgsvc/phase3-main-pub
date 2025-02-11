
SON_TO_FATHER = {}
for line in open(config["reads"], 'r'):
	if line.startswith('#'):
		continue
	fields = line.strip().split()
	sex = fields[4]
	if (sex == '1') and (fields[2] != '0' ):
		SON_TO_FATHER[fields[1]] = fields[2]

print(SON_TO_FATHER)


rule extract_genotyped_sample:
	input:
		"{results}/population-typing/{callset}/{version}/{coverage}/merged-vcfs/filtered/all-samples_bi_all_{filter}.vcf.gz"
	output:
		temp("{results}/evaluation-chrY/{callset}/{version}/{coverage}/tmp/{sample}_bi_all_{filter}.vcf.gz")
	conda:
		"../envs/genotyping.yml"
	wildcard_constraints:
		filter = "unfiltered|strict|lenient"
	resources:
		mem_total_mb = 40000
	shell:
		"""
		bcftools view --samples {wildcards.sample} -O z -o {output}
		tabix -p vcf {output} 
		"""


rule extract_genotyped_variant_ids:
	input:
		vcf = "{results}/evaluation-chrY/{callset}/{version}/{coverage}/tmp/{sample}_bi_all_{filter}.vcf.gz",
		regions = regions_to_bed
	output:
		temp("{results}/evaluation-chrY/{callset}/{version}/{coverage}/tmp/{sample}_{filter}_{vartype}_{regions}.tsv")
	wildcard_constraints:
		filter = "unfiltered|strict|lenient"
	shell:
		"""
		zcat {input} | workflow/scripts/extract-varianttype.py {wildcards.vartype} | bedtools intersect -header -a - -b {input.regions} -u -f 0.5 | python3 workflow/scripts/get_ids.py > {output}
		"""


rule genotype_concordance_father_son:
	input:
		father = lambda wildcards: expand("{{results}}/evaluation-chrY/{{callset}}/{{version}}/{sample}/{{coverage}}/tmp/{sample}_bi_all_{{filter}}.vcf.gz", sample = SON_TO_FATHER[wildcards.son]),
		son = "{results}/evaluation-chrY/{callset}/{version}/{sample}/{coverage}/tmp/{sample}_bi_all_{filter}.vcf.gz",
		regions = regions_to_bed,
		typed_ids = "{results}/evaluation-chrY/{callset}/{version}/{coverage}/tmp/{sample}_{filter}_{vartype}_{regions}.tsv" 
	output:
		tmp_vcf1 = temp("{results}/evaluation-chrY/{callset}/{version}/{sample}/{coverage}/concordance/{regions}_{vartype}_base.vcf"),
		tmp_vcf2 = temp("{results}/evaluation-chrY/{callset}/{version}/{sample}{coverage}/concordance/{regions}_{vartype}_call.vcf"),
		summary = "{results}/evaluation-chrY/{callset}/{version}/{sample}/{coverage}/concordance/{regions}_{vartype}/summary.txt"
	conda:
		"../envs/genotyping.yml"
	log:
		"{results}/evaluation-chrY/{callset}/{version}/{sample}/{coverage}/concordance/{regions}_{vartype}/summary.log"
	resources:
		mem_total_mb = 40000,
		runtime_hrs = 0,
		runtime_min = 40
	shell:
		"""
		bedtools intersect -header -a {input.father} -b {input.regions} -u -f 0.5 | bgzip > {output.tmp_vcf1}
		bedtools intersect -header -a {input.son} -b {input.regions} -u -f 0.5 | bgzip > {output.tmp_vcf2}
		python3 workflow/scripts/genotype-evaluation.py {output.tmp_vcf1} {output.tmp_vcf2} {input.typed_ids} --qual 0 2> {log} 1> {output.summary}
		"""
	
	
rule collect_concordance_results:
	input:
		lambda wildcards: expand("{results}/evaluation-chrY/{{callset}}/{{version}}/{sample}/{{coverage}}/concordance/{{regions}}_{{vartype}}/summary.txt", sample=SON_TO_FATHER.keys())
	output:
		 "{results}/evaluation-chrY/{callset}/{version}/plots/{coverage}/concordance_{callset}-{version}-{coverage}_{regions}_{vartype}.tsv"
	params:
		samples = ','.join([c for c in SON_TO_FATHER.keys()]),
		outfile = "{results}/evaluation-chrY/{callset}/{version}/plots/{coverage}/{metric}_{callset}-{version}-{coverage}_{regions}_{vartype}",
		folder = "{results}/evaluation-chrY/{callset}/{version}"
	shell:
		"""
		python3 workflow/scripts/collect-results.py {wildcards.metric} {wildcards.coverage} {params.samples} {wildcards.regions} -variants {wildcards.vartype} -folder {params.folder} -outfile {params.outfile}
		"""
