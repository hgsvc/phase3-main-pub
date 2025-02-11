
# stores paths to reads
reads_leave_one_out = {}

for line in open(READS, 'r'):
	if line.startswith('#'):
		continue
	fields = line.strip().split()
	sample_name = fields[1]
	read_path = fields[7]
	reads_leave_one_out[sample_name] = read_path

allowed_variants = ['snp', 'indels', 'large-deletion', 'large-insertion', 'large-complex']
callsets_leave_one_out = [s for s in CALLSETS.keys()]
coverages_leave_one_out = ['full'] + [c for c in DOWNSAMPLING]
versions_leave_one_out = [v for v  in PANGENIE.keys()] + [v for v in PANGENIE_MODULES.keys()]


################################################################
######   prepare input panel and ground truth genotypes  #######
################################################################

# remove positions that are ".|." in the left out sample. These cannot be used for evaluation, as the true genotype is unknown
# for now, also remove CHM13, because PanGenie cannot handle haploid samples
# remove sample from the panel
rule remove_missing:
	input:
		lambda wildcards: CALLSETS[wildcards.callset]['multi'] if wildcards.representation == 'multi' else CALLSETS[wildcards.callset]['bi']
	output:
		temp("{results}/leave-one-out/{callset}/preprocessed-vcfs/{sample}_{callset}_{representation}_no-missing.vcf")
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_total_mb=20000
	priority: 1
	wildcard_constraints:
		representation = "bi|multi"
	shell:
		"zcat {input} | python3 workflow/scripts/remove-missing.py {wildcards.sample} > {output}"


rule prepare_panel:
	input:
		lambda wildcards: CALLSETS[wildcards.callset]['multi']
	output:
		temp("{results}/leave-one-out/{callset}/input-panel/panel-{sample}_{callset}.vcf")
	conda:
		"../envs/genotyping.yml"
	priority: 1
	log:
		"{results}/leave-one-out/{callset}/input-panel/panel-{sample}_{callset}.log"
	resources:
		mem_total_mb=20000
	shell:
		"bcftools view --samples ^{wildcards.sample} {input} | bcftools view --min-ac 1 2> {log} 1> {output}"



# extract ground truth genotypes for sample
rule prepare_truth:
	input:
		"{results}/leave-one-out/{callset}/preprocessed-vcfs/{sample}_{callset}_bi_no-missing.vcf.gz"
	output:
		temp("{results}/leave-one-out/{callset}/truth/truth-{sample}_{callset}.vcf")
	conda:
		"../envs/genotyping.yml"
	priority: 1
	resources:
		mem_total_mb=20000
	log:
		"{results}/leave-one-out/{callset}/truth/truth-{sample}_{callset}.log"
	shell:
		"bcftools view --samples {wildcards.sample} {input} 2> {log} 1> {output}"



rule compress_vcf:
	input:
		"{results}/leave-one-out/{filename}.vcf"
	output:
		vcf = "{results}/leave-one-out/{filename}.vcf.gz",
		tbi = "{results}/leave-one-out/{filename}.vcf.gz.tbi"
	priority: 1
	shell:
		"""
		bgzip -c {input} > {output.vcf}
		tabix -p vcf {output.vcf}
		"""


########################################################
##################    run PanGenie    ##################
########################################################


# run pangenie
rule pangenie:
	input:
		reads = lambda wildcards: reads_leave_one_out[wildcards.sample] if wildcards.coverage == 'full' else "{results}/downsampling/{callset}/{coverage}/{sample}_{coverage}.fa.gz",
		fasta = lambda wildcards: CALLSETS[wildcards.callset]['reference'],
		vcf="{results}/leave-one-out/{callset}/input-panel/panel-{sample}_{callset}.vcf"
	output:
		genotyping = temp("{results}/leave-one-out/{callset}/{version}/{sample}/{coverage}/temp/pangenie-{sample}_genotyping.vcf")
	log:
		"{results}/leave-one-out/{callset}/{version}/{sample}/{coverage}/pangenie-{sample}.log"
	threads: 24
	resources:
		mem_total_mb=190000,
		runtime_hrs=7,
		runtime_min=1
	priority: 1
	params:
		out_prefix="{results}/leave-one-out/{callset}/{version}/{sample}/{coverage}/temp/pangenie-{sample}",
		pangenie = lambda wildcards: PANGENIE[wildcards.version]
	wildcard_constraints:
		version = "|".join([k for k in PANGENIE.keys()] + ['^' + k for k in PANGENIE_MODULES])
	shell:
		"""
		(/usr/bin/time -v {params.pangenie} -i <(zcat {input.reads}) -v {input.vcf} -r /hilbert{input.fasta} -o {params.out_prefix} -s {wildcards.sample} -j {threads} -t {threads} -g ) &> {log}
		"""


# run pangenie in the modularized way (> v2.1.1)
rule pangenie_modules:
	input:
		reads = lambda wildcards: reads_leave_one_out[wildcards.sample] if wildcards.coverage == 'full' else "{results}/downsampling/{callset}/{coverage}/{sample}_{coverage}.fa.gz",
		fasta = lambda wildcards: CALLSETS[wildcards.callset]['reference'],
		vcf="{results}/leave-one-out/{callset}/input-panel/panel-{sample}_{callset}.vcf"
	output:
		genotyping = temp("{results}/leave-one-out/{callset}/{version}/{sample}/{coverage}/temp/pangenie-{sample}_genotyping.vcf"),
		index = temp(directory("{results}/leave-one-out/{callset}/{version}/{sample}/{coverage}/temp/index/"))
	log:
		index = "{results}/leave-one-out/{callset}/{version}/{sample}/{coverage}/pangenie-{sample}_index.log",
		genotype = "{results}/leave-one-out/{callset}/{version}/{sample}/{coverage}/pangenie-{sample}.log"
	threads: 24
	resources:
		mem_total_mb=190000,
		runtime_hrs=15,
		runtime_min=1
	priority: 1
	params:
		out_prefix="{results}/leave-one-out/{callset}/{version}/{sample}/{coverage}/temp/pangenie-{sample}",
		index_prefix="{results}/leave-one-out/{callset}/{version}/{sample}/{coverage}/temp/index/pangenie-{sample}",
		pangenie = lambda wildcards: PANGENIE_MODULES[wildcards.version].split('PanGenie')[0] + " PanGenie",
		pangenie_params = lambda wildcards: PANGENIE_MODULES[wildcards.version].split('PanGenie')[-1]
	wildcard_constraints:
		version = "|".join([k for k in PANGENIE_MODULES.keys()] + ['^' + k for k in PANGENIE])
	shell:
		"""
		mkdir {output.index}
		(/usr/bin/time -v {params.pangenie}-index -v {input.vcf} -r /hilbert{input.fasta} -o {params.index_prefix} -t {threads} ) &> {log.index}
		(/usr/bin/time -v {params.pangenie} {params.pangenie_params} -f {params.index_prefix} -i <(gunzip -c {input.reads}) -o {params.out_prefix} -j {threads} -t {threads} -s {wildcards.sample} ) &> {log.genotype}
		"""


########################################################
##################    Evaluation      ##################
########################################################


rule alleles_per_bubble:
	input:
		lambda wildcards: CALLSETS[wildcards.callset]['multi']
	output:
		plot = "{results}/leave-one-out/{callset}/alleles-per-bubble.pdf",
		bed = "{results}/leave-one-out/{callset}/complex-bubbles.bed"
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_total_mb=20000,
		runtime_hrs=1
	shell:
		"zcat {input} | python3 workflow/scripts/variant-statistics.py {output.plot} 1 > {output.bed}"


# prepare beds for biallelic and complex graph regions
rule prepare_beds:
	input:
		bed = "{results}/leave-one-out/{callset}/complex-bubbles.bed",
		fai = lambda wildcards: CALLSETS[wildcards.callset]['reference'] + '.fai'
	output:
		bed = "{results}/leave-one-out/{callset}/biallelic-bubbles.bed",
		tmp = temp("{results}/leave-one-out/{callset}/biallelic-bubbles.fai"),
		bed_tmp = temp("{results}/leave-one-out/{callset}/biallelic-bubbles.tmp")
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		sort -k1,1d -k 2,2n -k 3,3n {input.fai} > {output.tmp}
		sort -k1,1d -k 2,2n -k 3,3n {input.bed} > {output.bed_tmp}
		bedtools complement -i {output.bed_tmp} -g {output.tmp} > {output.bed}
		"""


# convert genotyped VCF to biallelic representation
rule convert_genotypes_to_biallelic:
	input:
		vcf = "{results}/leave-one-out/{callset}/{version}/{sample}/{coverage}/temp/pangenie-{sample}_genotyping.vcf.gz",
		biallelic = lambda wildcards: CALLSETS[wildcards.callset]['bi']
	output:
		temp("{results}/leave-one-out/{callset}/{version}/{sample}/{coverage}/pangenie-{sample}_genotyping-biallelic.vcf")
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_total_mb=30000
	priority: 1
	shell:
		"zcat {input.vcf} | python3 workflow/scripts/convert-to-biallelic.py {input.biallelic} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n \"}}' > {output}"


# determine untypable ids based on unmerged callsets (i.e. independent of graph construction)
# this does not account for IDs that possibly went missing during construction of the graph (i.e. the multi-allelic input VCF)
rule untypable_ids:
	input:
		lambda wildcards: CALLSETS[wildcards.callset]['bi']
	output:
		lists = "{results}/leave-one-out/{callset}/untypable-ids/{sample}-untypable.tsv",
		summary = temp("{results}/leave-one-out/{callset}/untypable-ids-{sample}.tsv")
	params:
		out = "{results}/leave-one-out/{callset}/untypable-ids"
	conda:
		"../envs/genotyping.yml"
	shell:
		"zcat {input} | python3 workflow/scripts/untypable-ids-single.py {params.out} {wildcards.sample} > {output.summary}"


# determine IDs untypable because they have been filtered out during graph construction.
# Whenever a graph is constructed by merging callset into graph, some IDs might go missing (due to conflicts or ".").
# This step catches such cases as well. Since graph VCF might contain uncovered IDs, some untypable IDs might not get caught,
# that is why files are combined with other untypables (above).
rule generally_untypable_ids:
	input:
		biallelic = lambda wildcards: CALLSETS[wildcards.callset]['bi'],
		multiallelic = "{results}/leave-one-out/{callset}/input-panel/panel-{sample}_{callset}.vcf",
		untypables = "{results}/leave-one-out/{callset}/untypable-ids-{sample}.tsv"
	output:
		"{results}/leave-one-out/{callset}/untypable-ids/{sample}-untypable-all.tsv"
	resources:
		mem_total_mb = 50000
	log:
		"{results}/leave-one-out/{callset}/untypable-ids/{sample}-untypable-all.log"
	shell:
		"python3 workflow/scripts/untypable-ids-general.py {input.biallelic} {input.multiallelic} 2> {log} | cat - {input.untypables} | sort | uniq > {output}" 



# determine untypable IDs
rule remove_untypable:
	input:
		vcf = "{results}/leave-one-out/{callset}/{path}{sample}{other}.vcf",
		ids = "{results}/leave-one-out/{callset}/untypable-ids/{sample}-untypable-all.tsv"
	output:
		vcf = "{results}/leave-one-out/{callset}/{path}{sample}{other}-typable-{vartype}.vcf.gz",
		tbi = "{results}/leave-one-out/{callset}/{path}{sample}{other}-typable-{vartype}.vcf.gz.tbi"
	wildcard_constraints:
		callset = "|".join([v for v in CALLSETS.keys()]),
		sample = "|".join([s for s in reads_leave_one_out.keys()]),
		vartype = "|".join(allowed_variants)
	resources:
		mem_total_mb = 20000,
		runtime_hrs = 1
	priority: 1
	shell:
		"""
		cat {input.vcf} | python3 workflow/scripts/skip-untypable.py {input.ids} | python3 workflow/scripts/extract-varianttype.py {wildcards.vartype} | bgzip -c > {output.vcf}
		tabix -p vcf {output.vcf}
		"""

rule rtg_format:
	input:
		lambda wildcards: CALLSETS[wildcards.callset]['reference']
	output:
		directory("{results}/leave-one-out/{callset}/SDF")
	resources:
		mem_total_mb=20000
	priority: 1
	shell:
		'rtg format -o {output} {input}'


def region_to_bed(wildcards):
	if wildcards.regions == "biallelic":
		return "{results}/leave-one-out/{callset}/biallelic-bubbles.bed".format(results=wildcards.results, callset=wildcards.callset)
	if wildcards.regions == "multiallelic":
		return "{results}/leave-one-out/{callset}/complex-bubbles.bed".format(results=wildcards.results, callset=wildcards.callset)
	if wildcards.regions in CALLSETS[wildcards.callset]["regions"]:
		return CALLSETS[wildcards.callset]["regions"][wildcards.regions]
	assert(False)


# precision-recall
rule vcfeval:
	input:
		callset = "{results}/leave-one-out/{callset}/{version}/{sample}/{coverage}/pangenie-{sample}_genotyping-biallelic-typable-{vartype}.vcf.gz",
		callset_tbi = "{results}/leave-one-out/{callset}/{version}/{sample}/{coverage}/pangenie-{sample}_genotyping-biallelic-typable-{vartype}.vcf.gz.tbi",
		baseline = "{results}/leave-one-out/{callset}/truth/truth-{sample}_{callset}-typable-{vartype}.vcf.gz",
		baseline_tbi = "{results}/leave-one-out/{callset}/truth/truth-{sample}_{callset}-typable-{vartype}.vcf.gz.tbi",
		regions = region_to_bed,
		sdf = "{results}/leave-one-out/{callset}/SDF"
	output:
		summary = "{results}/leave-one-out/{callset}/{version}/{sample}/{coverage}/precision-recall-typable/{regions}_{vartype}/summary.txt"
	conda:
		"../envs/genotyping.yml"
	priority: 1
	wildcard_constraints:
		sample = "|".join([s for s in reads_leave_one_out.keys()]),
		vartype = "|".join(allowed_variants),
		coverage = "|".join(coverages_leave_one_out)
	params:
		tmp = "{results}/leave-one-out/{callset}/{version}/{sample}/{coverage}/precision-recall-typable/{regions}_{vartype}_tmp",
		outname = "{results}/leave-one-out/{callset}/{version}/{sample}/{coverage}/precision-recall-typable/{regions}_{vartype}",
		which = "--all-records"
	resources:
		mem_total_mb = 30000,
		runtime_hrs = 1,
		runtime_min = 40
	shell:
		"""
		rtg vcfeval -b {input.baseline} -c {input.callset} -t {input.sdf} -o {params.tmp} --ref-overlap --evaluation-regions {input.regions} {params.which} --Xmax-length 30000 > {output.summary}.tmp
		mv {params.tmp}/* {params.outname}/
		mv {output.summary}.tmp {output.summary}
		rm -r {params.tmp}
		"""


# determine the variants that went into re-typing per category
rule collect_typed_variants:
	input:
		callset = "{results}/leave-one-out/{callset}/preprocessed-vcfs/{sample}_{callset}_bi_no-missing.vcf.gz",
		regions= region_to_bed,
		ids="{results}/leave-one-out/{callset}/untypable-ids/{sample}-untypable-all.tsv"
	output:
		"{results}/leave-one-out/{callset}/genotyped-ids/{sample}_{regions}_{vartype}.tsv"
	conda:
		"../envs/genotyping.yml"
	wildcard_constraints:
		sample = "|".join([s for s in reads_leave_one_out.keys()]),
		vartype = "|".join(allowed_variants)
	resources:
		mem_total_mb=50000
	priority: 1
	shell:
		"zcat {input.callset} | python3 workflow/scripts/skip-untypable.py {input.ids} | python3 workflow/scripts/extract-varianttype.py {wildcards.vartype} | bedtools intersect -header -a - -b {input.regions} -u -f 0.5 | python3 workflow/scripts/get_ids.py > {output}"


# compute concordances
rule genotype_concordances:
	input:
		callset = "{results}/leave-one-out/{callset}/{version}/{sample}/{coverage}/pangenie-{sample}_genotyping-biallelic-typable-{vartype}.vcf.gz",
		callset_tbi = "{results}/leave-one-out/{callset}/{version}/{sample}/{coverage}/pangenie-{sample}_genotyping-biallelic-typable-{vartype}.vcf.gz.tbi",
		baseline = "{results}/leave-one-out/{callset}/truth/truth-{sample}_{callset}-typable-{vartype}.vcf.gz",
		baseline_tbi = "{results}/leave-one-out/{callset}/truth/truth-{sample}_{callset}-typable-{vartype}.vcf.gz.tbi",
		regions = region_to_bed,
		typed_ids = "{results}/leave-one-out/{callset}/genotyped-ids/{sample}_{regions}_{vartype}.tsv"
	output:
		tmp_vcf1 = temp("{results}/leave-one-out/{callset}/{version}/{sample}/{coverage}/concordance/{regions}_{vartype}_base.vcf"),
		tmp_vcf2 = temp("{results}/leave-one-out/{callset}/{version}/{sample}/{coverage}/concordance/{regions}_{vartype}_call.vcf"),
		summary = "{results}/leave-one-out/{callset}/{version}/{sample}/{coverage}/concordance/{regions}_{vartype}/summary.txt"
	conda:
		"../envs/genotyping.yml"
	wildcard_constraints:
		sample = "|".join([s for s in reads_leave_one_out.keys()]),
		vartype = "|".join(allowed_variants),
		coverage = "|".join(coverages_leave_one_out)
	log:
		"{results}/leave-one-out/{callset}/{version}/{sample}/{coverage}/concordance/{regions}_{vartype}/summary.log"
	resources:
		mem_total_mb = 40000,
		runtime_hrs = 0,
		runtime_min = 40
	priority: 1
	shell:
		"""
		bedtools intersect -header -a {input.baseline} -b {input.regions} -u -f 0.5 | bgzip > {output.tmp_vcf1}
		bedtools intersect -header -a {input.callset} -b {input.regions} -u -f 0.5 | bgzip > {output.tmp_vcf2}
		python3 workflow/scripts/genotype-evaluation.py {output.tmp_vcf1} {output.tmp_vcf2} {input.typed_ids} --qual 0 2> {log} 1> {output.summary}
		"""


########################################################
##################     Plotting      ###################
########################################################


# collect results across all samples
rule collect_results:
	input:
		lambda wildcards: expand("{{results}}/leave-one-out/{{callset}}/{{version}}/{sample}/{{coverage}}/{{metric}}/{{regions}}_{{vartype}}/summary.txt", sample = CALLSETS[wildcards.callset]['leave_one_out_samples'])
	output:
		"{results}/leave-one-out/{callset}/{version}/plots/{coverage}/{metric}_{callset}-{version}-{coverage}_{regions}_{vartype}.tsv"
	params:
		samples = lambda wildcards: ','.join([c for c in CALLSETS[wildcards.callset]['leave_one_out_samples']]),
		outfile = "{results}/leave-one-out/{callset}/{version}/plots/{coverage}/{metric}_{callset}-{version}-{coverage}_{regions}_{vartype}",
		folder = "{results}/leave-one-out/{callset}/{version}"
	priority: 1
	shell:
		"python3 workflow/scripts/collect-results.py {wildcards.metric} {wildcards.coverage} {params.samples} {wildcards.regions} -variants {wildcards.vartype} -folder {params.folder} -outfile {params.outfile}"



# plot results across different regions for a version
rule plotting_regions:
	input:
		lambda wildcards: expand("{{results}}/leave-one-out/{{callset}}/{{version}}/plots/{{coverage}}/{m}_{{callset}}-{{version}}-{{coverage}}_{regions}_{vartype}.tsv", m = wildcards.metric if wildcards.metric != 'untyped' else 'concordance', regions=['biallelic', 'multiallelic'] + [r for r in CALLSETS[wildcards.callset]["regions"].keys()], vartype=CALLSETS[wildcards.callset]['variants'])
	output:
		"{results}/leave-one-out/{callset}/plots/comparison-regions/{metric}/{metric}_{coverage}_{version}.pdf"
	wildcard_constraints:
		metric="concordance|precision-recall-typable|untyped",
		coverage = "|".join(coverages_leave_one_out)
	priority: 1
	conda:
		"../envs/genotyping.yml"
	params:
		sources = lambda wildcards: ' '.join([wildcards.callset + '-' + wildcards.version + '-' + wildcards.coverage + '_' + r for r in ['biallelic', 'multiallelic'] + [r for r in CALLSETS[wildcards.callset]["regions"].keys()]])
	shell:
		"python3 workflow/scripts/plot-results.py -files {input} -outname {output} -sources {params.sources} -metric {wildcards.metric}"





# plot results of different subsampling runs
rule plotting_versions:
	input:
		lambda wildcards: expand("{{results}}/leave-one-out/{{callset}}/{version}/plots/{{coverage}}/{m}_{{callset}}-{version}-{{coverage}}_{{regions}}_{vartype}.tsv", m = wildcards.metric if wildcards.metric != 'untyped' else 'concordance', version=versions_leave_one_out, vartype=CALLSETS[wildcards.callset]['variants'])
	output:
		"{results}/leave-one-out/{callset}/plots/comparison-versions/{metric}/{metric}_{coverage}_{regions}.pdf"
	wildcard_constraints:
		metric="concordance|precision-recall-typable|untyped",
		coverage = "|".join(coverages_leave_one_out)
	priority: 1
	conda:
		"../envs/genotyping.yml"
	params:
		sources = lambda wildcards: ' '.join([wildcards.callset + '-' + v + '-' + wildcards.coverage + '_' + wildcards.regions for v in versions_leave_one_out])
	shell:
		"python3 workflow/scripts/plot-results.py -files {input} -outname {output} -sources {params.sources} -metric {wildcards.metric}"



# plot results of different subsampling runs, comparing concordance and typed variants per sample
rule plotting_versions_conc_vs_untyped:
	input:
		lambda wildcards: expand("{{results}}/leave-one-out/{{callset}}/{version}/plots/{{coverage}}/concordance_{{callset}}-{version}-{{coverage}}_{{regions}}_{vartype}.tsv", version=versions_leave_one_out, vartype=CALLSETS[wildcards.callset]['variants'])
	output:
		"{results}/leave-one-out/{callset}/plots/comparison-versions/concordance-vs-untyped/concordance-vs-untyped_{coverage}_{regions}.pdf"
	wildcard_constraints:
		coverage = "|".join(coverages_leave_one_out)
	priority: 1
	conda:
		"../envs/genotyping.yml"
	params:
		sources = lambda wildcards: ' '.join([wildcards.callset + '-' + v + '-' + wildcards.coverage + '_' + wildcards.regions for v in versions_leave_one_out])
	shell:
		"python3 workflow/scripts/plot-results.py -files {input} -outname {output} -sources {params.sources} -metric concordance-vs-untyped"



# plot results of different coverages
rule plotting_coverages:
	input:
		lambda wildcards: expand("{{results}}/leave-one-out/{{callset}}/{{version}}/plots/{coverage}/{m}_{{callset}}-{{version}}-{coverage}_{{regions}}_{vartype}.tsv", m = wildcards.metric if wildcards.metric != 'untyped' else 'concordance', coverage=coverages_leave_one_out, vartype=CALLSETS[wildcards.callset]['variants'])
	output:
		"{results}/leave-one-out/{callset}/plots/comparison-coverages/{metric}/{metric}_{version}_{regions}.pdf"
	wildcard_constraints:
		metric="concordance|precision-recall-typable|untyped"
	priority: 1
	conda:
		"../envs/genotyping.yml"
	params:
		sources = lambda wildcards: ' '.join([wildcards.callset + '-' + wildcards.version + '-' + c + '_' + wildcards.regions for c in coverages_leave_one_out])
	shell:
		"python3 workflow/scripts/plot-results.py -files {input} -outname {output} -sources {params.sources} -metric {wildcards.metric}"



# plot resources (single core CPU time and max RSS) for different subsampling runs
rule plotting_resources:
	input:
		lambda wildcards: expand("{{results}}/leave-one-out/{{callset}}/{version}/{sample}/{{coverage}}/pangenie-{sample}.log", version = versions_leave_one_out, sample = CALLSETS[wildcards.callset]['leave_one_out_samples'])
	output:
		"{results}/leave-one-out/{callset}/plots/resources/resources_{callset}-{coverage}.pdf"
	conda:
		"../envs/genotyping.yml"
	params:
		outname = "{results}/leave-one-out/{callset}/plots/resources/resources_{callset}-{coverage}",
		samples	= lambda wildcards: " ".join(CALLSETS[wildcards.callset]['leave_one_out_samples']),
		versions = " ".join(versions_leave_one_out)
	shell:
		"python3 workflow/scripts/plot-resources.py -files {input} -outname {params.outname} -samples {params.samples} -sizes {params.versions}"
