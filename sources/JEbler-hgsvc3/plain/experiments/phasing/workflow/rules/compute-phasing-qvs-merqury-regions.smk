configfile: "config/config.yaml"

from random import shuffle

sample_to_reads = {}
for line in open(READS, 'r'):
	if line.startswith('#'):
		continue
	fields = line.strip().split()
	sample_to_reads[fields[1]] = fields[7]

chroms = []
for line in open(REFERENCE, 'r'):
	if line.startswith('>'):
		fields = line.strip().split()
		chrom = fields[0][1:]
		if not chrom in ["Y", "chrY", "M", "chrM"]:
			chroms.append(chrom)

chroms = {}
for callset in PHASED_VCFS:
	chroms[callset] = []
	for line in open(PHASED_VCFS[callset]["reference"]):
		if line.startswith('>'):
			fields = line.strip().split()
			chrom = fields[0][1:]
			if not chrom in ["Y", "chrY", "M", "chrM"] and (not "random" in chrom) and not ("chrUn" in chrom):
				chroms[callset].append(chrom)

samples = [s for s in QV_SAMPLES]
permuted_samples = [samples[i] for i in range(1, len(samples))] + [samples[0]]
randomized_samples = {samples[i] : permuted_samples[i] for i in range(len(samples))}

kmer_size = 21

qv_regions = set([])
for c in QV_REGIONS:
	for r in QV_REGIONS[c]:
		qv_regions.add(r)


rule compute_complement:
	"""
	This produces BED files defining all regions outside of the
	region of interest.
	"""
	input:
		bed = lambda wildcards: QV_REGIONS[wildcards.callset][wildcards.region],
		fai = lambda wildcards: PHASED_VCFS[wildcards.callset]["reference"] + ".fai"
	output:
		tmp_bed = temp("{results}/haplotypes/{callset}/{region}/bed/{region}_uncallable_tmp.bed"),
		tmp_fai = temp("{results}/haplotypes/{callset}/{region}/bed/{region}_genome.fai"),
		bed = "{results}/haplotypes/{callset}/{region}/bed/{region}_uncallable.bed"
	conda:
		"../envs/shapeit.yaml"
	shell:
		"""
		cat {input.bed} | sort -k1,1 -k2,2n > {output.tmp_bed}
		sort -k1,1 -k2,2n {input.fai} > {output.tmp_fai} 
		bedtools complement -i {output.tmp_bed} -g {output.tmp_fai} > {output.bed}
		"""


rule mask_genome:
	"""
	Set all bases outside of the desired regions to N.
	"""
	input:
		reference = lambda wildcards: PHASED_VCFS[wildcards.callset]["reference"],
		bed = "{results}/haplotypes/{callset}/{region}/bed/{region}_uncallable.bed"
	output:
		fasta = temp("{results}/haplotypes/{callset}/{region}/genome/genome_{region}_masked.fasta"),
		fasta_gz = temp("{results}/haplotypes/{callset}/{region}/genome/genome_{region}_masked.fasta.gz")
	log:
		"{results}/haplotypes/{callset}/{region}/genome/genome_{region}_masked.log"
	conda:
		"../envs/shapeit.yaml"
	resources:
		mem_total_mb = 10000
	shell:
		"""
		bedtools maskfasta -fi {input.reference} -fo {output.fasta} -bed {input.bed} &> {log}
		samtools faidx {output.fasta}
		bgzip -c {output.fasta} > {output.fasta_gz}
		"""


rule subset_vcf:
	"""
	Extract only QV samples from VCF to reduce file size.
	This helps speeding up next steps.
	"""
	input:
		lambda wildcards: PHASED_VCFS[wildcards.callset]["vcf"]
	output:
		temp("{results}/haplotypes/{callset}/samples.vcf.gz")
	log:
		"{results}/haplotypes/{callset}/samples.log"
	conda:
		"../envs/shapeit.yaml"
	threads:
		24
	resources:
		mem_total_mb = 20000,
		runtime_hrs = 4
	params:
		samples = ",".join(QV_SAMPLES)
	shell:
		"""
		bcftools view --threads {threads} -s {params.samples} {input} -O z -o {output} &> {log}
		tabix -p vcf {output}
		"""



rule compute_consensus_region:
	"""
	Insert all variants into the reference genome to produce haplotypes.
	"""
	input:
		vcf = "{results}/haplotypes/{callset}/samples.vcf.gz",
		reference = "{results}/haplotypes/{callset}/{region}/genome/genome_{region}_masked.fasta",
		bed = "{results}/haplotypes/{callset}/{region}/bed/{region}_uncallable.bed"
	output:
		fasta_gz = temp("{results}/haplotypes/{callset}/{region}/{callset}_{region}_{sample}_hap{haplotype}_consensus.fasta.gz"),
		tmp = temp("{results}/haplotypes/{callset}/{region}/{callset}_{region}_{sample}_hap{haplotype}_tmp.vcf.gz")
	log:
		"{results}/haplotypes/{callset}/{region}/{callset}_{region}_{sample}_hap{haplotype}.log"
	conda:
		"../envs/shapeit.yaml"
	wildcard_constraints:
		haplotype = "1|2",
		region = "|".join(qv_regions)
	resources:
		mem_total_mb = 20000,
		runtime_hrs = 4
	shell:
		"""
		bedtools subtract -header -A -a {input.vcf} -b {input.bed} | bgzip > {output.tmp}
		tabix -p vcf {output.tmp}		
		bcftools consensus --sample {wildcards.sample} --haplotype {wildcards.haplotype} -e 'ALT~\"<.*>\"' -f {input.reference} {output.tmp} 2> {log} | bgzip > {output.fasta_gz}
		"""


rule compute_consensus:
	"""
	Insert all variants into the reference genome to produce haplotypes.
	"""
	input:
		vcf = lambda wildcards: "{results}/haplotypes/{callset}/samples.vcf.gz" if wildcards.sample in QV_SAMPLES else PHASED_VCFS[wildcards.callset]["vcf"],
		reference = lambda wildcards: PHASED_VCFS[wildcards.callset]["reference"]
	output:
		fasta_gz = "{results}/haplotypes/{callset}/all/{callset}_all_{sample}_hap{haplotype}_consensus.fasta.gz"
	log:
		"{results}/haplotypes/{callset}/all/{callset}_all_{sample}_hap{haplotype}.log"
	conda:
		"../envs/shapeit.yaml"
	wildcard_constraints:
		haplotype = "1|2"
	resources:
		mem_total_mb = 20000,
		runtime_hrs = 4
	shell:
		"""
		bcftools consensus --sample {wildcards.sample} --haplotype {wildcards.haplotype} -e 'ALT~\"<.*>\"' -f {input.reference} {input.vcf} 2> {log} | bgzip > {output.fasta_gz}
		"""


rule extract_chromosomes_and_remove_masked:
	"""
	Exclude chrY, since it is not present
	in all samples.
	"""
	input:
		"{results}/haplotypes/{callset}/{region}/{callset}_{region}_{sample}_hap{haplotype}_consensus.fasta.gz"
	output:
		fasta = temp("{results}/haplotypes/{callset}/{region}/{callset}_{region}_{sample}_hap{haplotype}.fasta"),
		fasta_gz = temp("{results}/haplotypes/{callset}/{region}/{callset}_{region}_{sample}_hap{haplotype}.fasta.gz")
	log:
		"{results}/haplotypes/{callset}/{region}/{callset}_{region}_{sample}_hap{haplotype}.log"
	conda:
		"../envs/shapeit.yaml"
	wildcard_constraints:
		haplotype = "1|2"
	params:
		regions = lambda wildcards: " ".join(chroms[wildcards.callset]),
	shell:
		"""
		samtools faidx {input} {params.regions} | python3 workflow/scripts/split_fasta.py -o {output.fasta} &> {log}
		bgzip -c {output.fasta} > {output.fasta_gz}
		"""

rule meryl_counting:
	"""
	Counting kmers using meryl
	"""
	input:
		lambda wildcards: sample_to_reads[wildcards.sample]
	output:
		directory("{results}/read-counts-meryl/{sample}.meryl")
	threads:
		24
	resources:
		mem_total_mb = 32768,
		mem_total_gb = 32,
		time_hrs = 2
	conda:
		"../envs/merqury.yaml"
	log:
		"{results}/read-counts-meryl/{sample}.log"
	shell:
		"""
		meryl count k={kmer_size} memory={resources.mem_total_gb} threads={threads} {input} output {output} &> {log}
		"""



rule merqury_evaluate_assemblies:
	"""
	Compute assembly statistics using merqury.
	"""
	input:
		assembly = "{results}/haplotypes/{callset}/{region}/{callset}_{region}_{sample}_hap{haplotype}.fasta.gz",
		counts = "{results}/read-counts-meryl/{sample}.meryl"
	output:
		qv = "{results}/evaluation/assigned/{callset}/{region}/{sample}_hap{haplotype}/{callset}_{region}_{sample}_hap{haplotype}.qv",
		completeness = "{results}/evaluation/assigned/{callset}/{region}/{sample}_hap{haplotype}/{callset}_{region}_{sample}_hap{haplotype}.completeness.stats"
	threads:
		24
	resources:
		mem_total_mb = 32768,
		time_hrs = 2,
		mem_total_gb = 32
	wildcard_constraints:
		haplotype = "1|2"
	conda:
		"../envs/merqury.yaml"
	log:
		"{results}/evaluation/assigned/{callset}/{region}/{sample}_hap{haplotype}/{callset}_{region}_{sample}_hap{haplotype}.log"
	params:
		out_prefix = "{results}/evaluation/assigned/{callset}/{region}/{sample}_hap{haplotype}/{callset}_{region}_{sample}_hap{haplotype}"
	shell:
		"""
		./workflow/scripts/compute-qv.sh {input.counts} {input.assembly} {params.out_prefix} {threads} {resources.mem_total_gb} &> {log}
		./workflow/scripts/compute-completeness.sh {input.counts} {input.assembly} {params.out_prefix} {threads} {resources.mem_total_gb} 
		"""



rule merqury_evaluate_assemblies_randomized:
	"""
	Compute assembly statistics using merqury.
	"""
	input:
		assembly = "{results}/haplotypes/{callset}/{region}/{callset}_{region}_{sample}_hap{haplotype}.fasta.gz",
		counts = lambda wildcards: "{results}/read-counts-meryl/{s}.meryl".format(results=wildcards.results, s = randomized_samples[wildcards.sample])
	output:
		qv = "{results}/evaluation/randomized/{callset}/{region}/{sample}_hap{haplotype}/{callset}_{region}_{sample}_hap{haplotype}.qv",
		completeness = "{results}/evaluation/randomized/{callset}/{region}/{sample}_hap{haplotype}/{callset}_{region}_{sample}_hap{haplotype}.completeness.stats"
	threads:
		24
	resources:
		mem_total_mb = 32768,
		mem_total_gb = 32,
		time_hrs = 2
	conda:
		"../envs/merqury.yaml"
	log:
		"{results}/evaluation/randomized/{callset}/{region}/{sample}_hap{haplotype}/{callset}_{region}_{sample}_hap{haplotype}.log"
	params:
		out_prefix = "{results}/evaluation/randomized/{callset}/{region}/{sample}_hap{haplotype}/{callset}_{region}_{sample}_hap{haplotype}"
	shell:
		"""
		./workflow/scripts/compute-qv.sh {input.counts} {input.assembly} {params.out_prefix} {threads} {resources.mem_total_gb} &> {log}
		./workflow/scripts/compute-completeness.sh {input.counts} {input.assembly} {params.out_prefix} {threads} {resources.mem_total_gb} 
		"""


def merqury_plot_qv_input(wildcards):
	files = []
	for callset in PHASED_VCFS.keys():
		files.extend( expand("{results}/evaluation/{mode}/{callset}/{region}/{sample}_hap{haplotype}/{callset}_{region}_{sample}_hap{haplotype}.qv", results=wildcards.results, mode=["assigned", "randomized"], sample=QV_SAMPLES, callset=callset, haplotype=["1", "2"], region=["all"] + [r for r in QV_REGIONS[callset].keys()]) )
	return files



rule merqury_plot_qv:
	"""
	Plot merged + single QV values.
	"""
	input:
		computed_qvs = merqury_plot_qv_input,
		given_qvs = ASSEMBLY_QVS
	output:
		"{results}/evaluation/qv-values_merqury.pdf"
	log:
		"{results}/evaluation/qv-values_merqury.log"
	conda:
		"../envs/plotting.yaml"
	shell:
		"""
		ls {input.computed_qvs} | python3 workflow/scripts/plot-qv-merqury.py {output} -assembly {input.given_qvs} &> {log}
		"""



def merqury_plot_completeness_input(wildcards):
	files = []
	for callset in PHASED_VCFS.keys():
		files.extend( expand("{results}/evaluation/{mode}/{callset}/{region}/{sample}_hap{haplotype}/{callset}_{region}_{sample}_hap{haplotype}.completeness.stats", results=wildcards.results, mode=["assigned", "randomized"], sample=QV_SAMPLES, callset=callset, haplotype=["1", "2"], region=["all"] + [r for r in QV_REGIONS[callset].keys()]) )
	return files





rule merqury_plot_completeness:
	"""
	Plot merged + single completeness values.
	"""
	input:
		computed_qvs = merqury_plot_completeness_input
	output:
		"{results}/evaluation/completeness_merqury.pdf"
	log:
		"{results}/evaluation/completeness_merqury.log"
	conda:
		"../envs/plotting.yaml"
	shell:
		"""
		ls {input.computed_qvs} | python3 workflow/scripts/plot-completeness-merqury.py {output} &> {log}
		"""


rule compress_haplotypes:
	input:
		genomes = expand("{{results}}/haplotypes/{{callset}}/all/{{callset}}_all_{sample}_hap{haplotype}_consensus.fasta.gz", sample = CONSENSUS_SAMPLES, haplotype = ["1", "2"]),
		reference = lambda wildcards: PHASED_VCFS[wildcards.callset]["reference"]
	output:
		"{results}/haplotypes/{callset}/all/{callset}_all_consensus-haplotypes.agc"
	log:
		"{results}/haplotypes/{callset}/all/{callset}_all_consensus-haplotypes.log"
	conda:
		"../envs/agc.yaml"
	resources:
		mem_total_mb = 30000,
		runtime_hrs = 20
	threads: 32
	shell:
		"""
		agc create {input.reference} {input.genomes} -o {output} -t {threads} &> {log}
		"""


rule compress_haplotypes_qv:
	input:
		genomes = expand("{{results}}/haplotypes/{{callset}}/all/{{callset}}_all_{sample}_hap{haplotype}_consensus.fasta.gz", sample = QV_SAMPLES, haplotype = ["1", "2"]),
		reference = lambda wildcards: PHASED_VCFS[wildcards.callset]["reference"]
	output:
		"{results}/haplotypes/{callset}/all/{callset}_qv_consensus-haplotypes.agc"
	log:
		"{results}/haplotypes/{callset}/all/{callset}_qv_consensus-haplotypes.log"
	conda:
		"../envs/agc.yaml"
	resources:
		mem_total_mb = 30000,
		runtime_hrs = 20
	threads: 32
	shell:
		"""
		agc create {input.reference} {input.genomes} -o {output} -t {threads} &> {log}
		"""


rule collect_stats:
	input:
		expand("{{results}}/evaluation/assigned/{{callset}}/all/{sample}_hap{haplotype}/{{callset}}_all_{sample}_hap{haplotype}.qv", sample=QV_SAMPLES, haplotype=["1", "2"])
	output:
		"{results}/haplotypes/{callset}/all/{callset}_all_summary.tsv"
	shell:
		"""
		ls {input} | python3 workflow/scripts/collect-stats.py > {output}
		"""
