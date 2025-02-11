

rule extract_from_agc:
	input:
		AGC
	output:
		"{results}/data/{sample}_{haplotype}/{sample}_{haplotype}_consensus.fa.gz"
	conda:
		"../envs/agc.yaml"
	params:
		name = lambda wildcards: CONSENSUS_NAMES[wildcards.sample][wildcards.haplotype]
	resources:
		mem_mb = 20000
	shell:
		"""
		agc getset {input} {params.name} | bgzip > {output}
		samtools faidx {output}
		"""
