from collections import defaultdict

samples = defaultdict(lambda: dict())

for line in open(READS, 'r'):
	if line.startswith('#'):
		continue
	fields = line.split()
	sample = fields[1]
	reads = fields[7]
	samples[sample]['full'] = reads

	for coverage in DOWNSAMPLING:
		samples[sample][coverage] = "{results}/downsampling/" + coverage + "/" + sample + "_" + coverage + ".fa.gz"


regions = defaultdict(list)
chromosomes = defaultdict(list)
chromosome_lengths = defaultdict(lambda: defaultdict(int))

for callset in CALLSETS.keys():
	for line in open(CALLSETS[callset]['reference'] + '.fai', 'r'):
		fields = line.strip().split()
		chrom = fields[0]
		if not chrom in [str(i) for i in range(1,23)] + ['X', 'Y'] and not chrom in ['chr' + str(i) for i in range(1,23)] + ['chrX', 'chrY']:
			continue
		chromosomes[callset].append(chrom)
		chromosome_lengths[callset][chrom] = int(fields[1])

step_size = 1000000
for callset in CALLSETS.keys():
	for chromosome in sorted(chromosomes[callset]):
		pos = 1
		while pos < chromosome_lengths[callset][chromosome]:
			end_pos = min(pos+step_size-1, chromosome_lengths[callset][chromosome])
			regions[callset].append('{}:{}-{}'.format(chromosome, pos, end_pos))
			pos = pos + step_size

print('Running population genotyping on ' + str(len(samples)) + ' samples.')

##################################################
#####   genotype variants using PanGenie    ######
##################################################


# if there are variants with AF=0, remove them
rule prepare_population_panel:
	input:
		lambda wildcards: CALLSETS[wildcards.callset]['multi']
	output:
		"{results}/population-typing/{callset}/panel.vcf"
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_total_mb = 20000
	shell:
		"bcftools view --min-ac 1 {input} > {output}"


# genotype variants with pangenie
rule genotyping:
	input:
		reads=lambda wildcards: samples[wildcards.sample][wildcards.coverage],
		reference=lambda wildcards: CALLSETS[wildcards.callset]['reference'],
		panel = "{results}/population-typing/{callset}/panel.vcf"
	output:
		reads = temp("{results}/population-typing/{callset}/{version}/{coverage}/genotyping/reads-{sample}.fa"),
		path_segments = temp("{results}/population-typing/{callset}/{version}/{coverage}/genotyping/pangenie-{sample}_path_segments.fasta"),
		genotyping_vcf = temp("{results}/population-typing/{callset}/{version}/{coverage}/genotyping/pangenie-{sample}_genotyping.vcf"),
	log:
		"{results}/population-typing/{callset}/{version}/{coverage}/genotyping/pangenie-{sample}.log"
	wildcard_constraints:
		version = "|".join([k for k in PANGENIE.keys()] + ['^' + k for k in PANGENIE_MODULES])
	threads: 24
	resources:
		mem_total_mb = lambda wildcards: 180000 if "v100" in wildcards.version else 100000,
		runtime_hrs = 5,
		runtime_min = 1
	params:
		out_prefix = "{results}/population-typing/{callset}/{version}/{coverage}/genotyping/pangenie-{sample}",
		pangenie = lambda wildcards: PANGENIE[wildcards.version]
	benchmark:
		"{results}/population-typing/{callset}/{version}/{coverage}/genotyping/pangenie-{sample}_benchmark.txt"
	shell:
		"""
		gunzip -c {input.reads} > {output.reads}
		(/usr/bin/time -v {params.pangenie} -i {output.reads} -v {input.panel} -r /hilbert{input.reference} -o {params.out_prefix} -s {wildcards.sample} -j {threads} -t {threads} -g ) &> {log}
		"""


# genotype variants with pangenie in the modularized way - indexing step (to be done once)
rule genotyping_index:
	input:
		reference=lambda wildcards: CALLSETS[wildcards.callset]['reference'],
		panel = "{results}/population-typing/{callset}/panel.vcf"
	output:
		directory("{results}/population-typing/{callset}/{version}/{coverage}/genotyping/indexing/")
	log:
		"{results}/population-typing/{callset}/{version}/{coverage}/genotyping/indexing/pangenie-index.log"
	wildcard_constraints:
		version = "|".join([k for k in PANGENIE_MODULES.keys()] + ['^' + k for k in PANGENIE])
	threads: 24
	resources:
		mem_total_mb = 80000,
		runtime_hrs = 5,
		runtime_min = 1
	params:
		out_prefix = "{results}/population-typing/{callset}/{version}/{coverage}/genotyping/indexing/index",
		pangenie = lambda wildcards: PANGENIE_MODULES[wildcards.version].split('PanGenie')[0] + "PanGenie"
	benchmark:
		"{results}/population-typing/{callset}/{version}/{coverage}/genotyping/indexing/indexing_benchmark.txt"
	shell:
		"""
		(/usr/bin/time -v {params.pangenie}-index -v {input.panel} -r /hilbert{input.reference} -o {params.out_prefix} -t {threads} ) &> {log}
		"""


# genotype variants with pangenie in modularized way - genotyping step (to be done per sample)
rule genotyping_genotype:
	input:
		directory("{results}/population-typing/{callset}/{version}/{coverage}/genotyping/indexing/"),
		reads=lambda wildcards: samples[wildcards.sample][wildcards.coverage]
	output:
		genotyping_vcf = temp("{results}/population-typing/{callset}/{version}/{coverage}/genotyping/pangenie-{sample}_genotyping.vcf"),
	log:
		"{results}/population-typing/{callset}/{version}/{coverage}/genotyping/pangenie-{sample}.log"
	wildcard_constraints:
		 version = "|".join([k for k in PANGENIE_MODULES.keys()] + ['^' + k for k in PANGENIE])
	threads: 24
	resources:
		mem_total_mb = 75000,
		runtime_hrs = 5,
		runtime_min = 1
	params:
		out_prefix = "{results}/population-typing/{callset}/{version}/{coverage}/genotyping/pangenie-{sample}",
		in_prefix = "{results}/population-typing/{callset}/{version}/{coverage}/genotyping/indexing/index",
		pangenie = lambda wildcards: PANGENIE_MODULES[wildcards.version].split('PanGenie')[0] + " PanGenie",
		pangenie_params = lambda wildcards: PANGENIE_MODULES[wildcards.version].split('PanGenie')[-1]
	benchmark:
		"{results}/population-typing/{callset}/{version}/{coverage}/genotyping/pangenie-{sample}_benchmark.txt"
	shell:
		"""
		(/usr/bin/time -v {params.pangenie} {params.pangenie_params} -f {params.in_prefix} -i <(gunzip -c {input.reads}) -o {params.out_prefix} -j {threads} -t {threads} -s {wildcards.sample} ) &> {log}
		"""





# compress multiallelic file and create version without SNVs
rule postprocess_multi:
	input:
		"{results}/population-typing/{callset}/{version}/{coverage}/genotyping/pangenie-{sample}_genotyping.vcf"
	output:
		"{results}/population-typing/{callset}/{version}/{coverage}/genotyping/{sample}_genotyping_multi_all.vcf.gz"
	log:
		"{results}/population-typing/{callset}/{version}/{coverage}/genotyping/{sample}_genotyping_postprocess.log"
	resources:
		mem_total_mb=10000,
		runtime_hrs=0,
		runtime_min=50
	conda:
		"../envs/genotyping.yml"
	benchmark:
		"{results}/population-typing/{callset}/{version}/{coverage}/genotyping/{sample}_genotyping_postprocess_benchmark.txt"
	priority: 1
	shell:
		"(/usr/bin/time -v sh workflow/scripts/postprocess-multi.sh {input} {output} ) &> {log}"


# tabix vcf file
rule tabix_vcf:
	input:
		"{results}/population-typing/{filename}.vcf.gz"
	output:
		"{results}/population-typing/{filename}.vcf.gz.tbi"
	resources:
		mem_total_mb=20000,
		runtime_hrs=2,
		runtime_min=59
	conda:
		"../envs/genotyping.yml"
	shell:
		"tabix -p vcf {input}"
		

ruleorder: create_biallelic_vcf > tabix_vcf


# prepare list of multiallelic files to be merged
rule create_file_list_multi:
	input:
		vcfs=lambda wildcards: expand("{{results}}/population-typing/{{callset}}/{{version}}/{{coverage}}/genotyping/{sample}_genotyping_multi_all.vcf.gz", sample=sorted(samples.keys())),
#		tbi=lambda wildcards: expand("{{results}}/population-typing/{{callset}}/{{version}}/{{coverage}}/genotyping/{sample}_genotyping_multi_all.vcf.gz.tbi", sample=sorted(samples.keys()))
	output:
		"{results}/population-typing/{callset}/{version}/{coverage}/multi_all_filelist.tsv"
	run:
		f = open(output[0], 'w')
		for name in input.vcfs:
			print(name, file=f)
		f.close()


# merge single-sample vcfs into multi-sample vcf for a specific region
rule merge_vcfs_by_region_multi:
	input:
		"{results}/population-typing/{callset}/{version}/{coverage}/multi_all_filelist.tsv"
	output:
		"{results}/population-typing/{callset}/{version}/{coverage}/merged-vcfs/region-wise/pangenie_merged_multi_all_{region}.vcf.gz"
	log:
		"{results}/population-typing/{callset}/{version}/{coverage}/merged-vcfs/region-wise/pangenie_merged_multi_all_{region}.log"
	benchmark:
		"{results}/population-typing/{callset}/{version}/{coverage}/merged-vcfs/region-wise/pangenie_merged_multi_all_{region}_benchmark.txt"
	resources:
#		mem_total_mb = 50000,
		mem_total_mb = 800000,
		runtime_hrs = 2,
		runtime_min = 59
	wildcard_constraints:
		region = "chr[0-9A-Z]+:[0-9]+-[0-9]+"
	conda:
		"../envs/whatshap.yml" # NOTE: bcftools.yml has a problem using +fill-tags plugin
	shell:
		"(/usr/bin/time -v bcftools merge -r {wildcards.region} -l {input} | bcftools +fill-tags -Oz -o {output} -- -t AN,AC,AF) &> {log} "


# prepare list of biallelic files to be merged
rule create_file_list_bi:
	input:
		vcfs = lambda wildcards: expand("{{results}}/population-typing/{{callset}}/{{version}}/{{coverage}}/genotyping/{sample}_genotyping_bi_all.vcf.gz", sample=sorted(samples.keys())),
#		tbi = lambda wildcards: expand("{{results}}/population-typing/{{callset}}/{{version}}/{{coverage}}/genotyping/{sample}_genotyping_bi_all.vcf.gz.tbi", sample=sorted(samples.keys()))
	output:
		"{results}/population-typing/{callset}/{version}/{coverage}/bi_all_filelist.tsv"
	run:
		f = open(output[0], 'w')
		for name in input.vcfs:
			print(name, file=f)
		f.close()


# merge all biallelic vcfs into a multisample vcf file
rule merge_vcfs_by_region_bi:
	input:
		"{results}/population-typing/{callset}/{version}/{coverage}/bi_all_filelist.tsv"
	output:
		"{results}/population-typing/{callset}/{version}/{coverage}/merged-vcfs/region-wise/pangenie_merged_bi_all_{region}.vcf.gz"
	log:
		"{results}/population-typing/{callset}/{version}/{coverage}/merged-vcfs/region-wise/pangenie_merged_bi_all_{region}.log"
	resources:
		mem_total_mb = 500000,
#		mem_total_mb = 1000000,
#		runtime_hrs = 7,
		runtime_hrs = 2,
		runtime_min = 59
	wildcard_constraints:
		region = "chr[0-9A-Z]+:[0-9]+-[0-9]+"
	threads: 12
	benchmark:
		"{results}/population-typing/{callset}/{version}/{coverage}/merged-vcfs/region-wise/pangenie_merged_bi_all_{region}_benchmark.txt"
	conda:
		"../envs/whatshap.yml" # NOTE: bcftools.yml has a problem using +fill-tags plugin
	shell:
		"(/usr/bin/time -v bcftools merge -r {wildcards.region} -l {input} --threads {threads} | bcftools +fill-tags -Oz -o {output} -- -t AN,AC,AF) &> {log} "


# combine all region-specific vcfs into one whole genome one
rule concat_vcfs_filelist:
	input:
		vcfs = lambda wildcards: expand("{{results}}/population-typing/{{callset}}/{{version}}/{{coverage}}/merged-vcfs/region-wise/pangenie_merged_{{what}}_all_{region}.vcf.gz", region=regions[wildcards.callset]),
	#	tbi = lambda wildcards: expand("{{results}}/population-typing/{{callset}}/{{version}}/{{coverage}}/merged-vcfs/region-wise/pangenie_merged_{{what}}_all_{region}.vcf.gz.tbi", region=regions[wildcards.callset])
	output:
		"{results}/population-typing/{callset}/{version}/{coverage}/merged-vcfs/region-wise/concat_{what}_all_filelist.tsv"
	run:
		f = open(output[0], 'w')
		for name in input.vcfs:
			print(name, file=f)
		f.close()


# concat per-region vcfs
rule concat_vcfs:
	input:
		"{results}/population-typing/{callset}/{version}/{coverage}/merged-vcfs/region-wise/concat_{what}_all_filelist.tsv"
	output:
		vcf="{results}/population-typing/{callset}/{version}/{coverage}/merged-vcfs/whole-genome/all-samples_{what}_all.vcf.gz"
	log:
		"{results}/population-typing/{callset}/{version}/{coverage}/merged-vcfs/whole-genome/all-samples_{what}_all.log"
	benchmark:
		"{results}/population-typing/{callset}/{version}/{coverage}/merged-vcfs/whole-genome/all-samples_{what}_all_benchmark.txt"
	resources:
		mem_total_mb=50000,
		runtime_hrs=29,
		runtime_min=59
	conda:
		"../envs/genotyping.yml"
	shell:
		"bcftools concat -o {output.vcf} -O z -f {input} &> {log} "


# remove SNVs from vcf file
rule filter_snvs:
	input:
		vcf="{results}/population-typing/{callset}/{version}/{coverage}/merged-vcfs/whole-genome/all-samples_{what}_all.vcf.gz",
		tbi="{results}/population-typing/{callset}/{version}/{coverage}/merged-vcfs/whole-genome/all-samples_{what}_all.vcf.gz.tbi"
	output:
		"{results}/population-typing/{callset}/{version}/{coverage}/merged-vcfs/whole-genome/all-samples_{what}_nosnvs.vcf.gz"
	log:
		"{results}/population-typing/{callset}/{version}/{coverage}/merged-vcfs/whole-genome/all-samples_{what}_nosnvs.log"
	benchmark:
		"{results}/population-typing/{callset}/{version}/{coverage}/merged-vcfs/whole-genome/all-samples_{what}_nosnvs_benchmark.txt"
	resources:
		mem_total_mb=50000,
		runtime_hrs=27,
		runtime_min=59
	conda:
		"../envs/genotyping.yml"
	shell:
		"bcftools view --exclude-types snps -o {output} -O z {input.vcf} &> {log} "


# convert multi-allelic vcf to biallelic representation
rule create_biallelic_vcf:
	input:
		template = lambda wildcards: CALLSETS[wildcards.callset]['bi'],
		vcf="{results}/population-typing/{callset}/{version}/{coverage}/genotyping/{sample}_genotyping_multi_all.vcf.gz"
	output:
		vcf_bi_all="{results}/population-typing/{callset}/{version}/{coverage}/genotyping/{sample}_genotyping_bi_all.vcf.gz",
		tbi_bi_all="{results}/population-typing/{callset}/{version}/{coverage}/genotyping/{sample}_genotyping_bi_all.vcf.gz.tbi"
	log:
		"{results}/population-typing/{callset}/{version}/{coverage}/genotyping/{sample}_genotyping_postprocess_bi.log"
	benchmark:
		"{results}/population-typing/{callset}/{version}/{coverage}/genotyping/{sample}_genotyping_postprocess_bi_benchmark.txt"
	resources:
		mem_total_mb=30000,
		runtime_hrs=1,
		runtime_min=59
	conda:
		"../envs/genotyping.yml"
	priority: 1
	shell:
		"(/usr/bin/time -v sh workflow/scripts/postprocess-bi.sh {input.template} {input.vcf} {output.vcf_bi_all}) &> {log}"
