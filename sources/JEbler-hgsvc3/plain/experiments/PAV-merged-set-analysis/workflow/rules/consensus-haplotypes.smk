configfile: "config/config.yaml"

single_calls = {}
callable_regions = {}
for line in open(config["single_vcfs"], 'r'):
	fields = line.strip().split()
	single_calls[fields[0]] = fields[1]
	callable_regions[(fields[0], '1')] = fields[2]
	callable_regions[(fields[0], '2')] = fields[3]

reads = {}
for line in open(config["reads"], 'r'):
	fields = line.strip().split()
	reads[fields[0]] = fields[1]

read_samples = [line.split()[0] for line in open(config["reads"], 'r')]
all_samples = [line.split()[0] for line in open(config["single_vcfs"], 'r')]
samples = [s for s in all_samples if s in read_samples]

yak = config['yak']



######################################
# prepare BED files for given regions
######################################

rule prepare_bed_easy:
	"""
	Prepare BED files: (callable_regions - region) intersected with giab_easy regions
	"""
	input:
		callable_bed = lambda wildcards: callable_regions[(wildcards.sample, wildcards.hap)],
		regions = lambda wildcards: config['bed'][wildcards.region],
#		easy = config['giab_easy_regions']
	output:
		result = "results/{region}/bed/regions/{sample}_{hap}_{region}.bed.gz",
#		tmp = temp("results/{region}/bed/regions/tmp_{sample}_{hap}_{region}.bed")
	conda:
		"../env/bcftools.yml"
	resources:
		mem_total_mb=50000
	shell:
		"""
		bedtools subtract -a {input.callable_bed} -b {input.regions} | bgzip > {output.result}
		"""



rule compute_complement:
	"""
	This produces BED files defining all regions outside of the
	callable regions.
	"""
	input:
		bed = lambda wildcards: callable_regions[(wildcards.sample, wildcards.hap)] if wildcards.region == "all" else "results/{region}/bed/regions/{sample}_{hap}_{region}.bed.gz",
		fai = config["reference"] + ".fai",
	output:
		tmp_bed = temp("results/{region}/bed/{sample}_{hap}_uncallable_tmp.bed"),
		tmp_fai = temp("results/{region}/bed/{sample}_{hap}_genome.bed"),
		bed = "results/{region}/bed/{sample}_{hap}_uncallable.bed"
	conda:
		"../env/bcftools.yml"
	shell:
		"""
		zcat {input.bed} | sort -k1,1 -k2,2n > {output.tmp_bed}
		sort -k1,1 -k2,2n {input.fai} > {output.tmp_fai} 
		bedtools complement -i {output.tmp_bed} -g {output.tmp_fai} > {output.bed}
		"""


rule mask_genome:
	"""
	Set all bases outside of the callable regions to N.
	"""
	input:
		reference = config["reference"],
		bed = "results/{region}/bed/{sample}_{hap}_uncallable.bed"
	output:
		fasta = temp("results/{region}/genome/genome_{sample}_hap{hap}_masked.fasta"),
		fasta_gz = temp("results/{region}/genome/genome_{sample}_hap{hap}_masked.fasta.gz")
	log:
		"results/{region}/genome/{sample}_{hap}_genome.log"
	conda:
		"../env/bcftools.yml"
	resources:
		mem_total_mb = 10000
	shell:
		"""
		bedtools maskfasta -fi {input.reference} -fo {output.fasta} -bed {input.bed} &> {log}
		samtools faidx {output.fasta}
		bgzip -c {output.fasta} > {output.fasta_gz}
		"""


rule compute_consensus:
	"""
	Insert all variants into the masked reference genome to produce haplotypes.
	"""
	input:
		vcf = lambda wildcards: config["merged_set"] if wildcards.calls == "merged" else single_calls[wildcards.sample],
		reference = "results/{region}/genome/genome_{sample}_hap{haplotype}_masked.fasta",
		uncallable = "results/{region}/bed/{sample}_{haplotype}_uncallable.bed"
	output:
		fasta = temp("results/{region}/{calls}/{calls}_{sample}_hap{haplotype}_masked.fasta.gz"),
		tmp = temp("results/{region}/{calls}/{calls}_{sample}_hap{haplotype}.vcf.gz"),
		tmp_tbi = temp("results/{region}/{calls}/{calls}_{sample}_hap{haplotype}.vcf.gz.tbi")
	log:
		"results/{region}/{calls}/{calls}_{sample}_hap{haplotype}_masked.log"
	conda:
		"../env/bcftools.yml"
	wildcard_constraints:
		calls = "merged|single",
		haplotype = "1|2"
	resources:
		mem_total_mb = 20000,
		runtime_hrs = 4
	shell:
		"""
		bedtools subtract -header -A -a {input.vcf} -b {input.uncallable} | bgzip > {output.tmp}
		tabix -p vcf {output.tmp}
		bcftools consensus --sample {wildcards.sample} --haplotype {wildcards.haplotype} -f {input.reference} {output.tmp} 2> {log} | bgzip > {output.fasta} 
		"""


rule remove_masked_regions:
	"""
	Remove regions outside of the callable regions. These
	have been masked by Ns in the earlier steps.
	"""
	input:
		"results/{region}/{calls}/{calls}_{sample}_hap{haplotype}_masked.fasta.gz"
	output:
		fasta_gz = "results/{region}/{calls}/{calls}_{sample}_hap{haplotype}.fasta.gz",
		fasta = temp("results/{region}/{calls}/{calls}_{sample}_hap{haplotype}.fasta")
	log:
		"results/{region}/{calls}/{calls}_{sample}_hap{haplotype}.log"
	wildcard_constraints:
		calls = "merged|single|genome",
		haplotype = "1|2"
	shell:
		"""
		zcat {input} | python3 workflow/scripts/split_fasta.py -o {output.fasta} > {log}
		bgzip -c {output.fasta} > {output.fasta_gz}
		"""



rule yak_counting:
	"""
	Counting step needed for QV estimation.
	"""
	input:
		lambda wildcards: reads[wildcards.sample]
	output:
		"results/read-counts/{sample}.yak"
	wildcard_constraints:
	threads:
		24
	resources:
		mem_total_mb = 80000,
		runtime_hrs = 3
	shell:
		"""
		{yak} count -k31 -b37 -t{threads} -o {output} {input}
		"""


rule evaluate_assemblies:
	"""
	Compute QV estimates using yak.
	"""
	input:
		assembly="results/{region}/{calls}/{calls}_{sample}_hap{haplotype}.fasta.gz",
		counts="results/read-counts/{sample}.yak"
	output:
		"results/{region}/evaluation/{calls}_{sample}_hap{haplotype}.txt"
	wildcard_constraints:
		calls = "merged|single|genome",
		haplotype = "1|2"
	resources:
		mem_total_mb = 80000,
		runtime_hrs = 4
	threads:
		24
	shell:
		"""
		{yak} qv -t{threads} -p -K3.2g -l100k {input.counts} {input.assembly} > {output}
		"""
