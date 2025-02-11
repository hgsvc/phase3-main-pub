
rule extract_region:
	input:
		bed = lambda wildcards: REGIONS[wildcards.region][wildcards.sample][wildcards.haplotype],
		vcf = "{results}/pav_{sample}_{haplotype}/{sample}_{haplotype}.vcf.gz"
	output:
		temp("{results}/evaluation/vcf/{sample}_{haplotype}_{region}.vcf.gz")
	conda:
		"../envs/plotting.yaml"
	shell:
		"""
		bedtools subtract -header -a {input.vcf} -b {input.bed} | bgzip > {output}
		"""


rule evaluate_calls:
	input:
		vcf = lambda wildcards: "{results}/pav_{sample}_{haplotype}/{sample}_{haplotype}.vcf.gz" if wildcards.region == "all" else "{results}/evaluation/vcf/{sample}_{haplotype}_{region}.vcf.gz",
		hap1 = lambda wildcards: REFERENCES[wildcards.sample]['H1'],
		hap2 = lambda wildcards: REFERENCES[wildcards.sample]['H2']
	output:
		tsv = "{results}/evaluation/{region}/{sample}_{haplotype}.tsv",
		bed_all = "{results}/evaluation/{region}/{sample}_{haplotype}_all.bed",
		bed_err = "{results}/evaluation/{region}/{sample}_{haplotype}_err.bed",
		intervals = "{results}/evaluation/{region}/{sample}_{haplotype}_intervals.tsv"
	params:
		qv1 = lambda wildcards: QV[wildcards.sample]['H1'],
		qv2 = lambda wildcards: QV[wildcards.sample]['H2'],
		ref = "consensus" if USE_CONS else "assembly"
	shell:
		"""
		zcat {input.vcf} | python3 workflow/scripts/evaluate-assemblies-intervals.py --name {wildcards.sample}_{wildcards.haplotype} --qv-h1 {params.qv1} --qv-h2 {params.qv2} --hap1 {input.hap1} --hap2 {input.hap2} --errors {output.bed_err} --all {output.bed_all} --reference {params.ref} --intervals {output.intervals} > {output.tsv}
		"""


rule aggregate_results:
	input:
		expand("{{results}}/evaluation/{{region}}/{sample}_{haplotype}.tsv", sample = REFERENCES.keys(), haplotype = ["H1", "H2"])
	output:
		"{results}/evaluation/{region}/all-stats.tsv"
	shell:
		"""
		head -n 1 {input[0]} > {output}; tail -n +2 -q {input} >> {output}
		"""


rule plot_variants:
	input:
		"{results}/evaluation/{region}/{sample}_{haplotype}_{which}.bed"
	output:
		"{results}/evaluation/{region}/plots/{sample}_{haplotype}_{which}.pdf"
	wildcard_constraints:
		which = "all|err"
	conda:
		"../envs/plotting.yaml"
	params:
		name = "{sample}_{haplotype}_{which}"
	resources:
		mem_mb = 10000	
	shell:
		"""
		cat {input} | python3 workflow/scripts/plot-svlen.py {params.name} {output}
		"""

rule plot_score_hist:
	input:
		"{results}/evaluation/{region}/{sample}_{haplotype}_intervals.tsv"
	output:
		"{results}/evaluation/{region}/plots/{sample}_{haplotype}_scores.pdf"
	conda:
		"../envs/plotting.yaml"
	log:
		"{results}/evaluation/{region}/plots/{sample}_{haplotype}_scores.log"
	params:
		name = "{sample}_{haplotype}"
	shell:
		"""
		cat {input} | python3 workflow/scripts/plot-score-hist.py {output} {params.name} &> {log}
		"""
