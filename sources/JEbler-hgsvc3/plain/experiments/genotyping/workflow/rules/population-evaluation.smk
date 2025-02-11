import gzip

cohort_samples = []
for line in open(READS, 'r'):
	if line.startswith('#'):
		continue
	fields = line.strip().split() 
	cohort_samples.append(fields[1])


rule collect_samples:
	input:
		READS
	output:
		all="{results}/population-typing/{callset}/sample-index.tsv",
		related="{results}/population-typing/{callset}/sample-index-related.tsv",
		unrelated="{results}/population-typing/{callset}/sample-index-unrelated.tsv"
	run:
		with open(output.all, 'w') as all_samples, open(output.related, 'w') as related_samples, open(output.unrelated, 'w') as unrelated_samples, open(input[0], 'r') as infile:
			for line in infile:
				if line.startswith("#"):
					continue
				fields = line.strip().split()
				sample = fields[1]
				all_samples.write(sample + '\n')
				if (fields[2] == '0') and (fields[3] == '0'):
					# unrelated sample
					unrelated_samples.write(sample + '\n')
				else:
					# related sample
					related_samples.write(sample + '\n')



####################################################################
# extract necessary subsets of samples (e.g. all unrelated samples)
####################################################################

rule extract_samples:
	input:
		vcf="{results}/population-typing/{callset}/{version}/{coverage}/merged-vcfs/whole-genome/all-samples_bi_all.vcf.gz",
		samples="{results}/population-typing/{callset}/sample-index-unrelated.tsv"
	output:
		"{results}/population-typing/{callset}/{version}/{coverage}/merged-vcfs/whole-genome/unrelated-samples_bi_all.vcf.gz"
	resources:
		mem_total_mb = 50000,
		runtime_hrs = 10
	conda:
		"../envs/genotyping.yml"
	shell:
		"bcftools view --samples-file {input.samples} --force-samples {input.vcf} | bgzip -c > {output}"



###################################################
# compute mendelian consistency
###################################################

# count variants mendelian consistent in 0,1,2,...,nr_trios trios
rule check_consistent_trios:
	input:
		vcf = "{results}/population-typing/{callset}/{version}/{coverage}/merged-vcfs/whole-genome/all-samples_bi_all.vcf.gz",
		ped = READS,
		samples = "{results}/population-typing/{callset}/sample-index.tsv"
	output:
		variant_stats="{results}/population-typing/{callset}/{version}/{coverage}/evaluation/statistics/all/mendelian-statistics_bi_all.tsv",
		trio_stats="{results}/population-typing/{callset}/{version}/{coverage}/evaluation/statistics/all/trio-statistics_bi_all.tsv"
	log:
		"{results}/population-typing/{callset}/{version}/{coverage}/evaluation/statistics/all/mendelian-statistics_bi_all.log"
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_total_mb=300000,
		runtime_hrs=96,
		runtime_min=59
	shell:
		"python3 workflow/scripts/mendelian-consistency.py statistics -vcf {input.vcf} -ped {input.ped} -samples {input.samples} -table {output.variant_stats} -column-prefix pangenie > {output.trio_stats}"


###################################################
# compute allele frequency/genotype statistics
###################################################

rule compute_statistics:
	input:
		vcf = "{results}/population-typing/{callset}/{version}/{coverage}/merged-vcfs/whole-genome/unrelated-samples_bi_all.vcf.gz",
		vcf_all = "{results}/population-typing/{callset}/{version}/{coverage}/merged-vcfs/whole-genome/all-samples_bi_all.vcf.gz",
		panel = lambda wildcards: CALLSETS[wildcards.callset]['bi']
	output:
		"{results}/population-typing/{callset}/{version}/{coverage}/evaluation/statistics/all/genotyping-statistics_bi_all.tsv"
	conda:
		'../envs/genotyping.yml'
	resources:
		mem_total_mb=200000,
		runtime_hrs=160,
		runtime_min=59
	shell:
		"python3 workflow/scripts/collect-vcf-stats.py {input.panel} {input.vcf} {input.vcf_all} > {output}"



#########################################
# self-genotyping evaluation
#########################################


# genotyping concordance for each ID (over all samples)
rule genotype_concordance_variants:
	input:
		computed="{results}/population-typing/{callset}/{version}/{coverage}/merged-vcfs/whole-genome/all-samples_bi_all.vcf.gz",
		true = lambda wildcards: CALLSETS[wildcards.callset]['bi']
	output:
		"{results}/population-typing/{callset}/{version}/{coverage}/evaluation/statistics/all/self_bi_all_variant-stats.tsv"
	params:
		file_prefix="{results}/population-typing/{callset}/{version}/{coverage}/evaluation/statistics/all/self_bi_all",
		column_prefix="pangenie_self-genotyping",
		# restrict to panel samples for which reads are available
		samples = lambda wildcards: ','.join( list( set(PANEL_SAMPLES[wildcards.callset]) & set(cohort_samples) ) )
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_total_mb=3000000,
		runtime_hrs=40,
		runtime_min=59
	log:
		"{results}/population-typing/{callset}/{version}/{coverage}/evaluation/statistics/all/self_bi_all_variant-stats.log"
	shell:
		"python3 workflow/scripts/genotype-concordance-variant.py {input.true} {input.computed} {params.file_prefix} {params.samples} {params.column_prefix} &> {log}"


#################################################
# make a table containing var ID and bubble ID
#################################################

rule variant_id_to_bubble:
	input:
		lambda wildcards: CALLSETS[wildcards.callset]['multi']
	output:
		"{results}/population-typing/{callset}/{version}/{coverage}/evaluation/statistics/all/bubble-statistics_bi_all.tsv"
	shell:
		"zcat {input} | python3 workflow/scripts/id_to_bubble.py > {output}"


#################################################
# annotate variants by BED file
#################################################

rule annotate_variants:
	input:
		vcf = lambda wildcards: CALLSETS[wildcards.callset]['bi'],
		bed = lambda wildcards: CALLSETS[wildcards.callset]["regions"][wildcards.regions]
	output:
		"{results}/population-typing/{callset}/{version}/{coverage}/evaluation/statistics/all/annotations_bi_all_{regions}.tsv"
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_total_mb=200000,
		runtime_hrs=14,
		runtime_min=59
	shell:
		"bedtools annotate -i {input.vcf} -files {input.bed} | python3 workflow/scripts/annotate_repeats.py -names {wildcards.regions} -format tsv > {output}"


#################################################
# plot the results
#################################################


rule merge_table:
	input:
		"{results}/population-typing/{callset}/{version}/{coverage}/evaluation/statistics/all/genotyping-statistics_bi_all.tsv",
		"{results}/population-typing/{callset}/{version}/{coverage}/evaluation/statistics/all/mendelian-statistics_bi_all.tsv",
		"{results}/population-typing/{callset}/{version}/{coverage}/evaluation/statistics/all/self_bi_all_variant-stats.tsv",
		"{results}/population-typing/{callset}/{version}/{coverage}/evaluation/statistics/all/bubble-statistics_bi_all.tsv",
#		lambda wildcards: expand("{{results}}/population-typing/{{callset}}/{{version}}/{{coverage}}/evaluation/statistics/all/annotations_bi_all_{regions}.tsv", regions=CALLSETS[wildcards.callset]["regions"].keys())
	output:
		"{results}/population-typing/{callset}/{version}/{coverage}/evaluation/statistics/all/summary_bi_all.tsv"
	conda:
		"../envs/plotting.yml"
	resources:
		mem_total_mb=150000,
		runtime_hrs=5,
		runtime_min=1
	shell:
		"python3 workflow/scripts/merge-tables.py {input} {output}"


rule perform_regression:
	input:
		"{results}/population-typing/{callset}/{version}/{coverage}/evaluation/statistics/all/summary_bi_all.tsv"
	output:
		filters = "{results}/population-typing/{callset}/{version}/{coverage}/evaluation/statistics/plot_bi_all_filters.tsv",
		regression = "{results}/population-typing/{callset}/{version}/{coverage}/evaluation/statistics/plot_bi_all_regression.tsv"
	params:
		outprefix="{results}/population-typing/{callset}/{version}/{coverage}/evaluation/statistics/plot_bi_all",
		threshold= 10 if len(cohort_samples) < 1000 else 50,
#		regions = lambda wildcards: "-r " + " ".join(CALLSETS[wildcards.callset]["regions"].keys()) if CALLSETS[wildcards.callset]["regions"] else ""
	log:
		"{results}/population-typing/{callset}/{version}/{coverage}/evaluation/statistics/plot_bi_all_regression.log"
	conda:
		'../envs/plotting.yml'
	resources:
		mem_total_mb=200000,
		runtime_hrs=15,
		runtime_min=59
	shell:
		"python3 workflow/scripts/analysis-stepwise.py -t {input} -o {params.outprefix} -n {params.threshold} --regression-only &> {log}"


rule plot_statistics:
	input:
		"{results}/population-typing/{callset}/{version}/{coverage}/evaluation/statistics/plot_bi_all_regression.tsv"
	output:
		"{results}/population-typing/{callset}/{version}/{coverage}/evaluation/statistics/plot_bi_all_plot.log"
	params:
		outprefix="{results}/population-typing/{callset}/{version}/{coverage}/evaluation/statistics/plot_bi_all",
		threshold= 10 if len(cohort_samples) < 1000 else 50,
#		regions = lambda wildcards: "-r " + " ".join(CALLSETS[wildcards.callset]["regions"].keys()) if CALLSETS[wildcards.callset]["regions"] else ""
	conda:
		'../envs/plotting.yml'
	resources:
		mem_total_mb=200000,
		runtime_hrs=15,
		runtime_min=59
	shell:
		"python3 workflow/scripts/analysis-stepwise.py -t {input} -o {params.outprefix} -n {params.threshold} --plot-only &> {output}"






################################################################
# compute filtered callset VCFs
################################################################

# compute filtered callsets
rule filtered_callsets:
	input:
		vcf="{results}/population-typing/{callset}/{version}/{coverage}/merged-vcfs/whole-genome/{population}_bi_all.vcf.gz",
		filters="{results}/population-typing/{callset}/{version}/{coverage}/evaluation/statistics/plot_bi_all_filters.tsv",
		fai= lambda wildcards: CALLSETS[wildcards.callset]['reference'] + ".fai"
	output:
		tmp = temp("{results}/population-typing/{callset}/{version}/{coverage}/merged-vcfs/filtered/tmp-{population}_bi_all_{filter}.vcf.gz"),
		vcf = "{results}/population-typing/{callset}/{version}/{coverage}/merged-vcfs/filtered/{population}_bi_all_{filter}.vcf.gz"
	resources:
		mem_total_mb=20000,
		runtime_hrs=10,
		runtime_min=59
	wildcard_constraints:
		filter="unfiltered|lenient|strict",
		population="all-samples|unrelated-samples"
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		zcat {input.vcf} | python3 workflow/scripts/select_ids.py {input.filters} {wildcards.filter} | bgzip -c > {output.tmp} 
		bcftools reheader --fai {input.fai} {output.tmp} > {output.vcf}
		tabix -p vcf {output.vcf}
		"""
