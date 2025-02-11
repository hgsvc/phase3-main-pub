

rule get_mc_filtered_regions:
	"""
	Determine which regions were filtered out from MC VCF,
	mainly due to being covered by less than 80% of haplotypes.
	"""
	input:
		unfiltered = RAW_MC_VCF,
		filtered = PANEL_MULTI
	output:
		"{results}/mc-filtered-records.bed"
	resources:
		mem_total_mb = 50000
	shell:
		"""
		python3 workflow/scripts/find-filtered-out-records.py {input.unfiltered} {input.filtered} > {output}
		"""


rule analyze_unmatched_pav:
	"""
	Check how many unmatched PAV calls fall into the regions
	filtered out from MC VCF (rule above).
	"""
	input:
		table = "{results}/{combination}/{combination}_intersections_{which}.tsv",
		bed = "{results}/mc-filtered-records.bed"
	output:
		unmatched_bed = "{results}/{combination}/{combination}_unmatched-pav_{which}.bed",
		result = "{results}/{combination}/{combination}_unmatched-pav_mc-missed_{which}.bed"
	params:
		column = lambda wildcards: CALLSETS[wildcards.combination]["pav_name"]
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		cat {input.table} | python3 workflow/scripts/get-unmatched.py -c {params.column} --write-bed > {output.unmatched_bed}
		bedtools subtract -a {output.unmatched_bed} -b {input.bed} -A > {output.result}
		"""



rule analyze_unmatched_graph:
	"""
	Extract MC variants unmatched in PAV calls and keep those,
	that have AF=0 across the HGSVC3 samples. 
	"""
	input:
		table = "{results}/{combination}/{combination}_intersections_{which}.tsv",
	 	vcf = PANEL_BI,
		samples = SAMPLES
	output:
		hprc_only = "{results}/{combination}/{combination}_unmatched-graph_hprc-only_{which}.bed",   # not in PAV and not called in HGSVC3 samples
		hgsvc_called = "{results}/{combination}/{combination}_unmatched-graph_hgsvc-called_{which}.bed", # not in PAV and called in HGSVC3 samples
		unmatched_mc = "{results}/{combination}/{combination}_unmatched-graph_{which}.bed"   # not in PAV
	wildcard_constraints:
		which = "jana|pete"
	conda:
		"../envs/bcftools.yaml"
	resources:
		mem_total_mb = 20000,
		runtime_hrs = 1
	params:
		column = lambda wildcards: CALLSETS[wildcards.combination]["graph_name"]
	shell:
		"""
		cat {input.table} | python3 workflow/scripts/get-unmatched.py -c {params.column} --write-bed > {output.unmatched_mc}
		bcftools view -S {input.samples} {input.vcf} | bcftools view --max-af 0.0 | python3 workflow/scripts/filter-alleles.py {output.unmatched_mc}  > {output.hprc_only}
		bcftools view -S {input.samples} {input.vcf} | bcftools view --min-ac 1 | python3 workflow/scripts/filter-alleles.py {output.unmatched_mc}  > {output.hgsvc_called}
		"""

rule define_end_regions:
	input:
		FAI
	output:
		"{results}/reference_end_regions.bed"
	shell:
		"awk '{{print $1\"\\t\"5000000\"\\t\"$2-5000000 }}' {input} | grep -v chrM > {output}"


rule analyze_end_regions:
	input:
		hgsvc_called = "{results}/{combination}/{combination}_unmatched-graph_hgsvc-called_{which}.bed",
		pav_mc_missed = "{results}/{combination}/{combination}_unmatched-pav_mc-missed_{which}.bed",
		end_regions = "{results}/reference_end_regions.bed"
	output:
		hgsvc_called = "{results}/{combination}/{combination}_unmatched-graph_hgsvc-called_telo_{which}.bed",
		pav_mc_missed = "{results}/{combination}/{combination}_unmatched-pav_mc-missed_telo_{which}.bed"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bedtools subtract -a {input.hgsvc_called} -b {input.end_regions} -wa > {output.hgsvc_called}
		bedtools subtract -a {input.pav_mc_missed} -b {input.end_regions} -wa > {output.pav_mc_missed}
		"""
		


