configfile: "config/config.yaml"

reads = {}
for line in open(config["reads"], 'r'):
        fields = line.strip().split()
        reads[fields[0]] = fields[1]

read_samples = [line.split()[0] for line in open(config["reads"], 'r')]
all_samples = [line.split()[0] for line in open(config["single_vcfs"], 'r')]
samples = [s for s in all_samples if s in read_samples]


rule plot_result_with_assemblies:
	"""
	Plot merged + single QV values.
	"""
	input:
		computed_qvs = expand( "results/{{region}}/evaluation/{calls}_{sample}_hap{haplotype}.txt", calls = ["merged", "single", "genome"], sample = samples, haplotype = [1,2]),
		given_qvs = config['assembly_qvs']
	output:
		"results/{region}/evaluation/qv-values_with-assemblies.pdf"
	conda:
		"../env/bcftools.yml"
	shell:
		"""
		ls {input.computed_qvs} | python3 workflow/scripts/plot-qv.py {output} -assembly {input.given_qvs}
		"""



rule plot_result_no_assemblies:
	"""
	Plot merged + single QV values.
	"""
	input:
		computed_qvs = expand( "results/{{region}}/evaluation/{calls}_{sample}_hap{haplotype}.txt", calls = ["merged", "single", "genome"], sample = samples, haplotype = [1,2])
	output:
		"results/{region}/evaluation/qv-values.pdf"
	conda:
		"../env/bcftools.yml"
	shell:
		"""
		ls {input.computed_qvs} | python3 workflow/scripts/plot-qv.py {output}
		"""
