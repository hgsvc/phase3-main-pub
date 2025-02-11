rule add_tags:
	input:
		lambda wildcards: CALLSETS[wildcards.combination]["graph_vcf"]
	output:
		"{results}/{combination}/graph-tagged.vcf"
	conda:
		"../envs/bcftools.yaml"
	resources:
		mem_total_mb = 10000,
		runtime_hrs = 5
	shell:
		"zcat {input} | python3 workflow/scripts/set-pass.py | bcftools view -G | python3 workflow/scripts/add-svtags.py > {output}"


def intersect_vcfs_files(wildcards):
	o_files = []
	files = []
	names = []
	if wildcards.combination == "all":
		for combination in CALLSETS.keys():
			filename = "{results}/{combination}/graph-tagged.vcf".format(results=wildcards.results, combination=combination)
			if not CALLSETS[combination]["graph_vcf"] in o_files:
				files.append(filename)
				names.append(CALLSETS[combination]["graph_name"])
				o_files.append(CALLSETS[combination]["graph_vcf"])
		for combination in CALLSETS.keys():
			filename = CALLSETS[combination]["pav_vcf"]
			if not filename in files:
				files.append(filename)
				names.append(CALLSETS[combination]["pav_name"])
	else:
		files.append("{results}/{combination}/graph-tagged.vcf".format(results=wildcards.results, combination=wildcards.combination))
		files.append(CALLSETS[wildcards.combination]["pav_vcf"])
		names.append(CALLSETS[wildcards.combination]["graph_name"])
		names.append(CALLSETS[wildcards.combination]["pav_name"])
	return files, names


def intersect_vcfs_names(wildcards):
	names = [] 
	if wildcards.combination == "all":
		for callset in ["graph", "pav"]:
			for combination in CALLSETS.keys():
				name = CALLSETS[combination][callset + '_name']
				if not name in names:
					names.append(name)
	else:
		names.append(CALLSETS[wildcards.combination]["graph_name"])
		names.append(CALLSETS[wildcards.combination]["pav_name"])
	return names


rule intersect_vcfs:
	input:
		lambda wildcards: intersect_vcfs_files(wildcards)[0]
	output:
		tsv = "{results}/{combination}/{combination}_intersections_jana.tsv",
		vcf = "{results}/{combination}/{combination}_intersections_jana.vcf",
		pdf = "{results}/{combination}/{combination}_distances_jana.pdf"
	conda:
		"../envs/upsetplot.yml"
	log:
		"{results}/{combination}/{combination}_intersections_jana.log"
	params:
		names = lambda wildcards: intersect_vcfs_files(wildcards)[1],
		columns = lambda wildcards: ["in_" + n for n in intersect_vcfs_files(wildcards)[1]]
	resources:
		mem_total_mb = 50000,
		runtime_hrs = 6
	shell:
		"""
		python3 workflow/scripts/intersect_callsets.py intersect -c {input} -n {params.names} -t {output.tsv} -v {output.vcf} -p {output.pdf} --id-from-vcf &> {log}
		"""


rule plot_intersections_jana:
	input:
		"{results}/{combination}/{combination}_intersections_jana.tsv"
	output:
		"{results}/{combination}/{combination}_intersections_jana.pdf"
	conda:
		"../envs/upsetplot.yml"
	log:
		"{results}/{combination}/{combination}_intersections_plot_jana.log"
	params:
		columns = lambda wildcards: ["in_" + n for n in intersect_vcfs_files(wildcards)[1]]
	resources:
		mem_total_mb = 10000,
		runtime_hrs = 1
	shell:
		"""
		python3 workflow/scripts/plot-upset.py -t {input} -o {output} -n {params.columns} &> {log}
		"""

rule analyze_gene_interruptive:
	input:
		tsv = "{results}/{combination}/{combination}_intersections_jana.tsv",
		links = LINKS,
		annotations = ANNOTATIONS
	output:
		"{results}/{combination}/{combination}_overlap-gene-interruptive_jana.tsv"
	params:
		columns = lambda wildcards: ["ID_" + n for n in intersect_vcfs_files(wildcards)[1]]
	shell:
		"""
		python3 workflow/scripts/match-pav-annotations.py {input.annotations} {input.links} {input.tsv} {params.columns} > {output}
		"""
