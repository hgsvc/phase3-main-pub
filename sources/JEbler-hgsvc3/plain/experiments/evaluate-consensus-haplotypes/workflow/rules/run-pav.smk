

rule prepare_assembly_table:
	input:
		cons1 = lambda wildcards: ASSEMBLIES[wildcards.sample]['H1'],
		cons2 = lambda wildcards: ASSEMBLIES[wildcards.sample]['H2'],
		assembly = lambda wildcards: REFERENCES[wildcards.sample]['H2'] if wildcards.haplotype == 'H1' else REFERENCES[wildcards.sample]['H1'],
		reference = lambda wildcards: REFERENCES[wildcards.sample][wildcards.haplotype]
	output:
		tsv = "{results}/pav_{sample}_{haplotype}/assemblies.tsv",
		h1 = temp("{results}/pav_{sample}_{haplotype}/h1.fa.gz"),
		h2 = temp("{results}/pav_{sample}_{haplotype}/h2.fa.gz"),
		h3 = temp("{results}/pav_{sample}_{haplotype}/h3.fa.gz"),
		ref = temp("{results}/pav_{sample}_{haplotype}/ref.fa.gz"),
		fai = temp("{results}/pav_{sample}_{haplotype}/ref.fa.gz.fai")
	shell:
		"""
		cp {input.cons1} {output.h1}
		cp {input.cons2} {output.h2}
		cp {input.assembly} {output.h3}
		cp {input.reference} {output.ref}
		cp {input.reference}.fai {output.fai}
		echo \"NAME\tHAP_h1\tHAP_h2\tHAP_h3\" > {output.tsv}
		echo \"{wildcards.sample}_{wildcards.haplotype}\th1.fa.gz\th2.fa.gz\th3.fa.gz\" >> {output.tsv}
		"""


rule prepare_pav_contig:
	input:
		fasta = "{results}/pav_{sample}_{haplotype}/ref.fa.gz"
	output:
		json = "{results}/pav_{sample}_{haplotype}/config.json"
	run:
		with open(output.json, 'w') as outfile:
			outfile.write("{\n")
			outfile.write("\t\"reference\": \"ref.fa.gz\",\n")
			outfile.write("\t\"no_link_qry\": \"True\"\n")
			outfile.write("}\n")


rule run_pav:
	input:
		tsv = "{results}/pav_{sample}_{haplotype}/assemblies.tsv",
		pav_config = "{results}/pav_{sample}_{haplotype}/config.json",
		h1 = "{results}/pav_{sample}_{haplotype}/h1.fa.gz",
		h2 = "{results}/pav_{sample}_{haplotype}/h2.fa.gz",
		h3 = "{results}/pav_{sample}_{haplotype}/h3.fa.gz",
		ref = "{results}/pav_{sample}_{haplotype}/ref.fa.gz",
		fai = "{results}/pav_{sample}_{haplotype}/ref.fa.gz.fai"
	output:
		complete = "{results}/pav_{sample}_{haplotype}/run.complete",
		vcf = "{results}/pav_{sample}_{haplotype}/{sample}_{haplotype}.vcf.gz"
	threads:
		24
	params:
		wdir = "{results}/pav_{sample}_{haplotype}/"
	log:
		"{results}/pav_{sample}_{haplotype}/pav.log"
	benchmark:
		"{results}/pav_{sample}_{haplotype}/pav.benchmark"
	resources:
		mem_mb = 150000,
		walltime = "24:00:00"	
	singularity:
		"workflow/container/pav_latest.sif"
	shell:
		"""
		snakemake --verbose --jobs {threads} -d {params.wdir} -s /opt/pav/Snakefile --rerun-incomplete --restart-times 0 --notemp &> {log} && touch {output.complete}
		"""
