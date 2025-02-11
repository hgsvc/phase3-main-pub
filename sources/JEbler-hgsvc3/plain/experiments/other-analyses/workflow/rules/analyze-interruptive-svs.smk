

rule liftover_variants:
	"""
	Liftover variants from GRCh38 to CHM13.
	Also, adds 10kb left and right of intervals
	"""
	input: 
		bed = lambda wildcards: REGIONS[wildcards.region],
		chain = CHAIN
	output:	
		tmp = "{results}/{region}-tmp.bed",
		unmapped = "{results}/{region}-unmapped.bed",
		lifted = "{results}/{region}.bed"
	shell:
		"""
		zcat {input.bed} | python3 workflow/scripts/add-tag.py > {output.tmp}
		{liftover} {output.tmp} {input.chain} {output.lifted} {output.unmapped} -multiple -bedPlus=3
		"""


rule prepare_header:
	output: "{results}/header.txt"
	shell:
		"""
		echo "##INFO=<ID=INTERVAL,Number=.,Type=String,Description=\\"Interval overlap.\\">" > {output}
		"""


rule filter_svs:
	"""
	Filter out SNPs and indels from the VCF.
	"""
	input:
		vcf = lambda wildcards: PANEL_BI if wildcards.callset == "panel" else PHASED_VCF,
		samples = SAMPLES
	output:
		vcf = temp("{results}/{callset}/variants-svs.vcf.gz"),
		tbi = temp("{results}/{callset}/variants-svs.vcf.gz.tbi")
	conda:
		"../envs/bcftools.yaml"
	wildcard_constraints:
		callset = "panel|phased"
	resources:
		mem_total_mb = 50000,
		runtime_hrs = 2
	shell:
		"""
		bcftools view -S {input.samples} {input.vcf} --force-samples | python3 workflow/scripts/extract-varianttype.py large | bgzip > {output.vcf}
		tabix -p vcf {output.vcf}
		"""


rule annotate_vcf:
	"""
	Annotates a VCF with a BED.
	IDs are generated from the columns of the BED file.
	"""
	input:
		bed = "{results}/{callset}/{region}.bed",
		vcf = "{results}/{callset}/variants-svs.vcf.gz",
		header = "{results}/header.txt" 
	output:
		bed = "{results}/{callset}/{region}-tagged.bed",
		vcf = "{results}/{callset}/annotated_{region}.vcf.gz"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		cat {input.bed} | python3 workflow/scripts/add-window.py > {output.bed}
		bcftools annotate {input.vcf} -a {output.bed} -c CHROM,FROM,TO,INTERVAL -h {input.header} -O z -o {output.vcf}
		"""


rule analyze_intervals:
	"""
	Given the annotations, count the number of variants
	contained in the intervals, possibly with the same type?
	"""
	input:
		vcf = "{results}/{callset}/annotated_{region}.vcf.gz",
		bed = "{results}/{callset}/{region}-tagged.bed",
		mc_missed = "results/jana-callset-comparison/mc-filtered-records.bed",
		filters = FILTERS,
		samples = SAMPLES
	output:
		tsv = "{results}/{callset}/{callset}_counts_{region}.tsv",
		tmp = "{results}/{callset}/{callset}_counts_{region}.txt"
	log:
		"{results}/{callset}/{callset}_counts_{region}.log"
	resources:
		mem_total_mb = 20000
	conda:
		"../envs/bcftools.yaml"
	params:
		outprefix = "{results}/{callset}/{callset}"
	shell:
		"""
		bedtools intersect -a {input.bed} -b {input.mc_missed} -wa -u  | cut -f 4 > {output.tmp}
		zcat {input.vcf} | python3 workflow/scripts/count-variants-in-intervals.py {input.filters} {input.bed} {input.samples} {params.outprefix} {output.tmp} 2> {log} 1> {output}
		"""

