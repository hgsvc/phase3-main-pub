

rule count_vars_in_regions:
	input:
		vcf = lambda wildcards: "{results}/population-typing/{callset}/{version}/{coverage}/merged-vcfs/whole-genome/{population}_bi_all.vcf.gz" if wildcards.filter == "unfiltered" else "{results}/population-typing/{callset}/{version}/{coverage}/merged-vcfs/filtered/{population}_bi_all_{filter}.vcf.gz",
		region = lambda wildcards: CALLSETS[wildcards.callset]["regions"][wildcards.regions]
	output:
		"{results}/population-typing/{callset}/{version}/{coverage}/evaluation/variant-counts/{population}_{filter}_{region}.tsv"
	resources:
		mem_total_mb = 20000,
		runtime_hrs = 1,
		runtime_min = 59
	wildcard_constraints:
		filter = "unfiltered|strict|lenient",
		population = "all-samples|unrelated-samples"
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		bedtools intersect -header -a {input.vcf} -b {input.region} -u -f 0.5 | python3 workflow/scripts/count-varianttype.py {wildcards.region} > {output}
		"""
