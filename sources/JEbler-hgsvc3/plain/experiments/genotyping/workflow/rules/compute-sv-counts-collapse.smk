callsets = [c for c in GENOTYPED_SETS.keys()]

## Note: this pipeline considers only chromosomes 1-22 + chr X

# also determine the list of samples contained in each VCF
rule prepare_vcf:
	input:
		lambda wildcards: GENOTYPED_SETS[wildcards.source]['vcf']
	output:
		samples = "{results}/{source}/vcfs/{source}_all.tsv"
	wildcard_constraints:
		source = '|'.join(callsets)
	log:
		"{results}/{source}/vcfs/{source}_all.log"
	conda:
		"../envs/truvari.yml"
	shell:
		"""
		bcftools query -l {input} > {output.samples}
		"""

# find out intersection of samples contained in all VCFs
rule intersect_samples:
	input:
		samples = expand("{{results}}/{source}/vcfs/{source}_all.tsv", source=callsets)
	output:
		"{results}/samples-intersection.tsv"
	run:
		result = None
		for sample in input.samples:
			names = set([s.strip() for s in open(sample, 'r')])
			if result is None:
				result = names
			else:
				result = result & names
		with open(output[0], 'w') as outfile:
			for s in result:
				outfile.write(s + '\n')


# keep only intersection of samples (so that VCFs are comparable). Merge similar alleles with truvari.
# add AF,AN,AC tags to VCFs (to make sure all VCFs have them).
rule extract_intersection_samples_collapse:
	input:
		vcf= lambda wildcards: GENOTYPED_SETS[wildcards.source]['vcf'],
		samples = "{results}/samples-intersection.tsv",
		reference = lambda wildcards: GENOTYPED_SETS[wildcards.source]['reference'],
		fai = lambda wildcards: GENOTYPED_SETS[wildcards.source]['reference'] + ".fai"
	output:
		unmerged=temp("{results}/{source}/vcfs/tmp-{source}_intersection_full-unmerged.vcf.gz"),
		merged="{results}/{source}/vcfs/{source}_intersection_full_collapsed.vcf.gz"
	conda:
		"../envs/truvari.yml"
	resources:
		mem_total_mb=80000,
		runtime_hrs=23,
		runtime_min=59
	threads: 10
	shell:
		"""
		bcftools view --threads {threads} --samples-file {input.samples} {input.vcf} | bcftools reheader --fai {input.fai} | python3 workflow/scripts/extract-varianttype.py large | bgzip -c > {output.unmerged}
		tabix -p vcf {output.unmerged}
		truvari collapse -r 500 -p 0.95 -P 0.95 -s 50 -S 100000 -f {input.reference} -i {output.unmerged} | bcftools sort | bcftools view --min-ac 1 | bcftools +fill-tags -Oz -o {output.merged} -- -t AN,AC,AF
		tabix -p vcf {output.merged}
		"""


# keep only intersection of samples (so that VCFs are comparable)
# add AF,AN,AC tags to VCFs (to make sure all VCFs have them)
rule extract_intersection_samples:
	input:
		vcf= lambda wildcards: GENOTYPED_SETS[wildcards.source]['vcf'],
		samples = "{results}/samples-intersection.tsv"
	output:
		"{results}/{source}/vcfs/{source}_intersection_full_raw.vcf.gz"
	conda:
		"../envs/truvari.yml"
	wildcard_constraints:
		source="|".join(callsets)
	resources:
		mem_total_mb=80000,
		runtime_hrs=23,
		runtime_min=59
	shell:
		"""
		bcftools view --samples-file {input.samples} {input.vcf} | python3 workflow/scripts/extract-varianttype.py large | bcftools view --min-ac 1 | bcftools +fill-tags -Oz -o {output} -- -t AN,AC,AF
		tabix -p vcf {output}
		"""



###################################################################################################################################################################################################
###############################  create violin plots showing the number of SVs per sample in the different callsets (for different allele frequency cutoffs)  #####################################
###################################################################################################################################################################################################


#rule annotate_calls:
#	input:
#		vcf="{filename}.vcf.gz",
#		regions = [v for k,v in REGIONS.items()]
#	output:
#		temp("{filename}_annotated.txt.gz")
#	conda:
#		"../envs/genotyping.yml"
#	resources:
#		mem_total_mb=80000
#	params:
#		names= [k for k,v in REGIONS.items()]
#	shell:
#		"bedtools annotate -i {input.vcf} -files {input.regions} | python3 workflow/scripts/annotate_repeats.py -vcf {input.vcf} -names {params.names} -format vcf | bgzip > {output}"



rule plot_sv_counts_filtered:
	input:
		raw = lambda wildcards: expand("{{results}}/{source}/vcfs/{source}_intersection_full_raw.vcf.gz", source = [s for s in callsets if not GENOTYPED_SETS[s]['collapse']]),
		collapse = lambda wildcards: expand("{{results}}/{source}/vcfs/{source}_intersection_full_collapsed.vcf.gz", source = [s for s in callsets if GENOTYPED_SETS[s]['collapse']]),
		populations = POPULATIONS
	output:
		"{results}/sv-count-comparison_all.pdf"
	log:
		"{results}/sv-count-comparison_all.log"
	conda:
		"../envs/plotting.yml"
	params:
		names_raw = ' '.join([s for s in callsets if not GENOTYPED_SETS[s]['collapse']]),
		names_collapsed = ' '.join([s for s in callsets if GENOTYPED_SETS[s]['collapse']]),
		outname = "{results}/sv-count-comparison_all",
	shell:
		"python3 workflow/scripts/plot-sv-counts.py -vcfs {input.collapse} {input.raw} -names {params.names_collapsed} {params.names_raw} -o {params.outname} -pop {input.populations} &> {log}"



rule plot_sv_counts_filtered_het:
	input:
		raw = lambda wildcards: expand("{{results}}/{source}/vcfs/{source}_intersection_full_raw.vcf.gz", source = [s for s in callsets if not GENOTYPED_SETS[s]['collapse']]),
		collapse = lambda wildcards: expand("{{results}}/{source}/vcfs/{source}_intersection_full_collapsed.vcf.gz", source = [s for s in callsets if GENOTYPED_SETS[s]['collapse']]),
		populations = POPULATIONS
	output:
		"{results}/sv-count-comparison_het.pdf"
	log:
		"{results}/sv-count-comparison_het.log"
	conda:
		"../envs/plotting.yml"
	params:
		names_raw = ' '.join([s for s in callsets if not GENOTYPED_SETS[s]['collapse']]),
		names_collapsed = ' '.join([s for s in callsets if GENOTYPED_SETS[s]['collapse']]),
		outname = "{results}/sv-count-comparison_het",
	shell:
		"python3 workflow/scripts/plot-sv-counts.py -vcfs {input.collapse} {input.raw} -names {params.names_collapsed} {params.names_raw} -o {params.outname} -pop {input.populations} --het &> {log}"





#rule plot_sv_counts_filtered_annotated:
#	input:
#		raw = lambda wildcards: expand("{{results}}/{source}/vcfs/{source}_intersection_full_raw_annotated.txt.gz", source = [s for s in callsets if not GENOTYPED_SETS[s]['collapse']]),
#		collapse = lambda wildcards: expand("{{results}}/{source}/vcfs/{source}_intersection_full_collapsed_annotated.txt.gz", source = [s for s in callsets if GENOTYPED_SETS[s]['collapse']]),
#		populations = POPULATIONS
#	output:
#		"{results}/sv-count-comparison_annotated.pdf"
#	log:
#		"{results}/sv-count-comparison_annotated.log"
#	conda:
#		"../envs/plotting.yml"
#	params:
#		names_raw = ' '.join([s for s in callsets if not GENOTYPED_SETS[s]['collapse']]),
#		names_collapsed = ' '.join([s for s in callsets if GENOTYPED_SETS[s]['collapse']]),
#		outname = "{results}/sv-count-comparison_annotated",
#		annotations =  [k for k,v in REGIONS.items()]
#	shell:
#		"python3 workflow/scripts/plot-sv-counts.py -vcfs {input.collapse} {input.raw} -names {params.names_collapsed} {params.names_raw} -o {params.outname} -pop {input.populations} --annotations {params.annotations} &> {log}"





######################################################################################################################
############################# plot length distribution for the different callsets ####################################
######################################################################################################################


rule plot_length_distribution:
	input:
		raw = lambda wildcards: expand("{{results}}/{source}/vcfs/{source}_intersection_full_raw.vcf.gz", source = [s for s in callsets if not GENOTYPED_SETS[s]['collapse']]),
		collapse = lambda wildcards: expand("{{results}}/{source}/vcfs/{source}_intersection_full_collapsed.vcf.gz", source = [s for s in callsets if GENOTYPED_SETS[s]['collapse']]),
		populations = POPULATIONS
	output:
		"{results}/length-distribution.pdf"
	log:
		"{results}/length-distribution.log"
	conda:
		"../envs/plotting.yml"
	params:
		names_raw = ' '.join([s for s in callsets if not GENOTYPED_SETS[s]['collapse']]),
		names_collapsed = ' '.join([s for s in callsets if GENOTYPED_SETS[s]['collapse']]),
	shell:
		"python3 workflow/scripts/plot-variant-length.py --callsets {input.collapse} {input.raw} --names {params.names_collapsed} {params.names_raw} -a 0.05 -o {output} &> {log}"

