import sys
import argparse
from collections import defaultdict
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import gzip
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


chromosomes = ['chr' + str(i) for i in range(1,23)] + ['chrX']

def parse_vcf(filename, callsetname, outfile, length, sample_to_pop, annotations, hap):
	counts = defaultdict(lambda: 0)
	n_ann = len(annotations)

	samples = []
	nr_vars = 0
	for line in gzip.open(filename, 'rt'):
		fields = line.split()
		if line.startswith('##'):
			continue
		if line.startswith('#'):
			if n_ann > 0:
				samples = fields[9:-n_ann]
			else:
				samples = fields[9:]
			continue
		# only keep SVs (>= 50bp)
		ref_allele = fields[3]
		alt_allele = fields[4]
		# make sure the VCF is biallelic
		assert len(alt_allele.split(',')) == 1

		if not fields[0] in chromosomes:
			continue

		# only consider SVs (alleles >= length), unless it is a symbolic VCF
		if (len(ref_allele) < length) and (len(alt_allele) < length) and (not alt_allele.startswith('<')):
			continue

		nr_vars += 1
		# find allele frequency
		info_field = { k.split('=')[0] : k.split('=')[1] for k in fields[7].split(';') if '=' in k }
		assert 'AF' in info_field
		print(line)
		allele_freq = float(info_field['AF'])

		af_category = ""
		if allele_freq < 0.01:
			af_category = '<1%'
		elif allele_freq <= 0.05:
			af_category = '1-5%'
		else:
			af_category = '>5%'

		# count SVs present in each sample
		gt_index = fields[8].split(':').index('GT')
		fields_to_consider = fields[9:-n_ann] if (n_ann > 0) else fields[9:]
		for genotype_info, sample in zip(fields_to_consider, samples):
			genotype = genotype_info.split(':')[gt_index]
			assert genotype in ['0', '1', '.', '0/0', '0/1', '1/0', '1/1', './1', '1/.', '0/.', './0', '0|0', '0|1', '1|0', '1|1', './.', '.', '.|.', '.|1', '1|.', '0|.', '.|0']

			if not hap:
				# count all present SVs
				nr_alts = 1 if '1' in genotype else 0
			else:
				# count only heterozygous SVs
				if genotype in ['1', '0/1', '1/0', '0|1', '1|0']:
					nr_alts = 1
				else:
					nr_alts = 0
			
			if nr_alts > 0:
				counts[(sample, af_category, 'all')] += nr_alts
				counts[(sample, 'all', 'all')] += nr_alts

				# indicates whether there was an overlap with any of the annotation regions
				any_annotation = False
	
				for i,annotation in enumerate(annotations):
					if float(fields[-n_ann+i]) >= 0.5:
						any_annotation = True
						counts[(sample, 'all', annotation)] += nr_alts
						counts[(sample, af_category, annotation)] += nr_alts
				if not any_annotation:
					counts[(sample, 'all', 'none')] += nr_alts
					counts[(sample, af_category, 'none')] += nr_alts


	for sample in samples:
		pop = sample_to_pop[sample]
		is_afr = 'AFR' if pop == 'AFR' else 'non-AFR'
		for annotation in ['all', 'none'] + annotations:
			for allele_freq in ['all', '<1%', '1-5%', '>5%']:
				var_count = counts[(sample, allele_freq, annotation)]
				outfile.write('\t'.join([sample, callsetname, str(var_count), allele_freq, annotation, sample_to_pop[sample], is_afr, callsetname + '-' + is_afr]) + '\n')

	sys.stderr.write('Processed callset ' + callsetname +  ' (' + str(nr_vars) + ' variants, ' + str(len(samples)) + ' samples).\n')

		

def run_plotting(vcfs, names, outname, length, populations, annotations, only_regions, hap):
	y_name = 'number of heterozygous SVs per sample' if hap else 'number of SVs per sample'
	sample_to_pop = parse_populations(populations)
	with open(outname + '.tsv', 'w') as outfile:
		# write header
		outfile.write('\t'.join(['sample_name', 'callset', y_name, 'allele frequency', 'region', 'superpopulation', 'population', 'callset-population']) + '\n')
		for v,n in zip(vcfs, names):
			parse_vcf(v, n, outfile, length, sample_to_pop, annotations, hap)

	# create violin plots
	df = pd.read_csv(outname + '.tsv', sep='\t')

	with PdfPages(outname +  '.pdf') as pdf:

		# (1) plot number of SVs for different allele frequency cutoffs in all regions

		df_plot = df[df['region'] == 'all']
		fig, ax = plt.subplots(figsize=(15,5))
		# Hide the right and top spines
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)

		# grouped violinplot
		sns.violinplot(ax=ax, data=df_plot, x='allele frequency', y=y_name, hue='callset', cut=0)
		pdf.savefig()
		plt.close()


		# (2) create same plot but split by AFR vs. non-AFR also

		fig, ax = plt.subplots(figsize=(15,5))
		# Hide the right and top spines
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)

		# grouped violinplot
		sns.violinplot(ax=ax, data=df_plot, x='allele frequency', y=y_name, hue='callset-population', cut=0)
		pdf.savefig()
		plt.close()


		# (3) same plot as before, but only <1% allele frequency
		df_plot_rare = df_plot[df_plot['allele frequency'] == '<1%']
		fig, ax = plt.subplots(figsize=(15,5))
		# Hide the right and top spines
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)

		# grouped violinplot
		sns.violinplot(ax=ax, data=df_plot_rare, x='allele frequency', y=y_name, hue='callset-population', cut=0)
		pdf.savefig()
		plt.close()


		# (4) plot number of SVs in different regions (all allele frequencies)

		df_plot = df[df['allele frequency'] == 'all']
		if only_regions:
			df_plot = df_plot[df_plot['region'] != 'all']
		else:
			df_plot = df_plot[df_plot['region'] != 'none']
		fig, ax = plt.subplots(figsize=(15,5))
		# Hide the right and top spines
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)

		# grouped violinplot
		fig.suptitle('all allele frequencies')
		sns.violinplot(ax=ax, data=df_plot, x='region', y=y_name, hue='callset', cut=0)
		pdf.savefig()
		plt.close()


		# (5) plot number of SVs in different regions for common variants (AF>5%)

		df_plot = df[df['allele frequency'] == '>5%']
		if only_regions:
			df_plot = df_plot[df_plot['region'] != 'all']
		else:
			df_plot = df_plot[df_plot['region'] != 'none']
		fig, ax = plt.subplots(figsize=(15,5))
		# Hide the right and top spines
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)

		# grouped violinplot
		fig.suptitle('common variants (AF > 5%)')
		sns.violinplot(ax=ax, data=df_plot, x='region', y=y_name, hue='callset', cut=0)
		pdf.savefig()
		plt.close()



		# (6) plot number of SVs in different regions for rare variants (AF<1%)

		df_plot = df[df['allele frequency'] == '<1%']
		if only_regions:
			df_plot = df_plot[df_plot['region'] != 'all']
		else:
			df_plot = df_plot[df_plot['region'] != 'none']
		fig, ax = plt.subplots(figsize=(15,5))
		# Hide the right and top spines
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)

		# grouped violinplot
		fig.suptitle('rare variants (AF < 1%)')
		sns.violinplot(ax=ax, data=df_plot, x='region', y=y_name, hue='callset', cut=0)
		pdf.savefig()
		plt.close()



		# (7) for each callset, create scatterplot showing AFR/non-AFR samples

		for callsetname in names:
			plt.figure()
			fig, ax = plt.subplots()
			fig.suptitle(callsetname)
			# Hide the right and top spines
			ax.spines['right'].set_visible(False)
			ax.spines['top'].set_visible(False)

			sns.scatterplot(ax=ax, data=df[(df.callset.isin([callsetname])) & (df['allele frequency'].isin(['all'])) & (df['region'].isin(['all']))], x='allele frequency', y=y_name, hue='population', alpha=0.6)
			pdf.savefig()
			plt.close()


		# (8) create scatterplot comparing the number of SVs per sample across pairs of callsets

		pop_to_color = {
			'AFR' : '#DB7D27',
			'AMR' : '#D72519',
			'EUR' : '#2D6F91',
			'EAS' : '#41A22F',
			'SAS' : '#782B8A'
		}
		
		# store (callset, sample) -> [SV count, population]
		callset_to_values = defaultdict(lambda: [None, None])
		# store callset -> [sample names]
		callset_to_samples = defaultdict(set)

		df_sub = df[(df['allele frequency'].isin(['all'])) & (df['region'].isin(['all']))]
		print('size_sub', len(df_sub))
		for sample, callset, counts, pop in zip(df_sub['sample_name'], df_sub['callset'], df_sub[y_name], df_sub['superpopulation']):
			callset_to_values[(callset, sample)] = [counts, pop]
			callset_to_samples[callset].add(sample)

		callset_names = [c for c in callset_to_samples.keys()]

		for callset1 in callset_names:
			for callset2 in callset_names:
				if callset1 == callset2:
					continue

				# find intersection of samples in these two callsets
				joint_samples = callset_to_samples[callset1] & callset_to_samples[callset2]

				# find out counts
				counts_c1 = []
				counts_c2 = []
				colors = []
				for sample in joint_samples:
					counts_c1.append(callset_to_values[(callset1, sample)][0])
					counts_c2.append(callset_to_values[(callset2, sample)][0])
					pop_c1 = callset_to_values[(callset1, sample)][1]
					pop_c2 = callset_to_values[(callset2, sample)][1]
					assert pop_c1 == pop_c2
					print(sample, callset1, callset2) 
					colors.append(pop_to_color[pop_c1])
				
				plt.figure()
				fig, ax = plt.subplots()
				fig.suptitle(y_name)
				# Hide the right and top spines
				ax.spines['right'].set_visible(False)
				ax.spines['top'].set_visible(False)

				scatter = ax.scatter(counts_c1, counts_c2, c=colors, alpha=0.3)
				ax.set_xlabel(callset1)
				ax.set_ylabel(callset2)


				# build legend
				custom_lines = [Line2D([0], [0], color=pop_to_color['AFR'], lw=4),
						Line2D([0], [0], color=pop_to_color['AMR'], lw=4),
						Line2D([0], [0], color=pop_to_color['EUR'], lw=4),
						Line2D([0], [0], color=pop_to_color['EAS'], lw=4),
						Line2D([0], [0], color=pop_to_color['SAS'], lw=4)]
				ax.legend(custom_lines, ['AFR', 'AMR', 'EUR', 'EAS', 'SAS'])
				pdf.savefig()
				plt.close()



def parse_populations(filename):
	sample_to_pop = defaultdict(lambda: 'unknown')
	for line in open(filename, 'r'):
		if line.startswith('#'):
			continue
		fields = line.strip().split()
		sample_name = fields[1]
		population = fields[6]
		sample_to_pop[sample_name] = population
	return sample_to_pop


if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='plot-sv-counts.py', description=__doc__)
	parser.add_argument('-vcfs', metavar='VCFs', nargs='+', required=True, help='Callsets to be considered.')
	parser.add_argument('-names', metavar='NAMES', nargs='+', required=True, help='Callset names (one for each callset)')
	parser.add_argument('-pop', metavar='POPULATION', required=True, help="tsv files specifying the population of each sample.")
	parser.add_argument('--annotations', metavar='ANNOTATIONS', nargs='+', default=[], help='Annotations contained in the input VCF file.')
	parser.add_argument('-o', metavar='OUTNAME', required= True, help='name prefix of the output files (<prefix>.pdf and <prefix>.tsv will be produced).')
	parser.add_argument('--length', metavar='LENGTH', help='Minimum allele length to consider', default=0, type=int)
	parser.add_argument('--only-regions', action='store_true', help='Only show anotation regions in count plots and skip all.', default=False)
	parser.add_argument('--het', action='store_true', help='Count number of heterozygous SVs only.', default=False)
	args = parser.parse_args()

	# make sure that the same number of callsets/names is provided
	assert len(args.vcfs) == len(args.names)
	run_plotting(args.vcfs, args.names, args.o, args.length, args.pop, args.annotations, args.only_regions, args.het)
