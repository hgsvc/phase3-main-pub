

rule preprocess_graph:
	input:
		lambda wildcards: CALLSETS[wildcards.combination]["graph_vcf"]
	output:
		"{results}/{combination}/{combination}_preprocessed_graph.vcf.gz"
	conda:
		"../envs/bcftools.yaml"
	resources:
		mem_total_mb = 50000,
		runtime_hrs = 10
	threads: 10
	shell:
		"""
		bcftools view --threads {threads} -G {input} -O z -o {output}
		tabix -p vcf {output}
		"""


rule count_overlaps_and_extract_unmatched:
	input:
		pav_del = lambda wildcards: CALLSETS[wildcards.combination]["pav_del"],
		pav_ins = lambda wildcards: CALLSETS[wildcards.combination]["pav_ins"],
		pav_vcf = lambda wildcards: CALLSETS[wildcards.combination]["pav_vcf"],
		mc_vcf = lambda wildcards: "{results}/{combination}/{combination}_preprocessed_graph.vcf.gz"
	output:
		tsv = "{results}/{combination}/{combination}_intersections_pete.tsv",
		unmatched = "{results}/{combination}/{combination}_unmatched-ids_pete.tsv"
	conda:
		"../envs/upsetplot.yml"
	log:
		"{results}/{combination}/{combination}_intersections_pete.log"
	params:
		pav_name = lambda wildcards: CALLSETS[wildcards.combination]["pav_name"],
		graph_name = lambda wildcards: CALLSETS[wildcards.combination]["graph_name"]
	resources:
		mem_total_mb = 50000,
		runtime_hrs = 5
	shell:
		"""
		zcat {input.pav_del} {input.pav_ins} | python3 workflow/scripts/analyze-overlaps.py {input.mc_vcf} {input.pav_vcf} {params.graph_name} {params.pav_name} {output.unmatched} | awk 'NR == 1; NR > 1 {{ print $0 | \"sort -k2,2V -k3,3n -k4,4n\" }}' 2> {log} 1> {output.tsv}
		"""


rule plot_intersections_pete:
	input:
		"{results}/{combination}/{combination}_intersections_pete.tsv"
	output:
		"{results}/{combination}/{combination}_intersections_pete.pdf"
	conda:
		"../envs/upsetplot.yml"
	params:
		pav_name = lambda wildcards: CALLSETS[wildcards.combination]["pav_name"],
		graph_name = lambda wildcards: CALLSETS[wildcards.combination]["graph_name"]
	log:
		"{results}/{combination}/{combination}_intersections_plot_pete.log"
	resources:
		mem_total_mb = 10000,
		runtime_hrs = 2
	shell:
		"""
		python3 workflow/scripts/plot-upset.py -t {input} -o {output} -n in_{params.graph_name} in_{params.pav_name} &> {log}
		"""

	


rule analyze_unmatched:
	input:
		unmatched = "{results}/{combination}/{combination}_unmatched-ids_pete.tsv",
	 	vcf = PANEL_BI,
		samples = SAMPLES
	output:
		"{results}/{combination}/{combination}_unmatched-filtered.vcf.gz"
	conda:
		"../envs/bcftools.yaml"
	resources:
		mem_total_mb = 20000,
		runtime_hrs = 1
	shell:
		"""
		bcftools view -S {input.samples} {input.vcf} | bcftools view --max-af 0.0 | python3 workflow/scripts/filter-alleles.py {input.unmatched} | bgzip > {output}
		"""
