#!/usr/bin/python

import argparse
import gzip
import sys
from collections import namedtuple
from collections import defaultdict

AlleleStats = namedtuple('AlleleStats','af ac an untyped')
GenotypeStats = namedtuple('GenotypeStats', 'heterozygosity het total')

def compute_allele_statistics(fields):
	"""
	Compute allele related statistics.
	"""
	an = 0
	ac = 0
	unknown = 0
	gt_index = fields[8].split(':').index('GT')
	for genotype in fields[9:]:
		gt = genotype.split(':')[gt_index]
		# assuming biallelic variants
		alleles = [gt[0], gt[-1]]
		for a in alleles:
			if a == '.':
				unknown += 1
				continue
			an += 1
			ac += int(a)
	af = ac / max(1.0, float(an))
	return AlleleStats(str(af), str(ac), str(an), str(unknown))	


def compute_genotype_statistics(fields, qualities=None):
	"""
	Compute genotype related statistics.
	"""
	counts = defaultdict(int)
	het_genotypes = 0
	total_genotypes = 0

	format_field = fields[8].split(':')

	gt_index = format_field.index('GT')
	gq_index = format_field.index('GQ') if qualities else None

	for genotype in fields[9:]:
		gt = genotype.split(':')[gt_index]
		# assuming biallelic variant
		alleles = [gt[0], gt[-1]]
		if not '.' in alleles:
			quality = int(genotype.split(':')[gq_index]) if qualities else None
			total_genotypes += 1
			if ('0' in alleles) and ('1' in alleles):
				het_genotypes += 1
			# read GQ
			if qualities is not None:
				for q in qualities:
					if quality >= q:
						counts[q] += 1
	genotype_stats = GenotypeStats( str(het_genotypes / max(1.0,float(total_genotypes))), str(het_genotypes), str(total_genotypes))
	return genotype_stats, counts


parser = argparse.ArgumentParser(prog='collect-vcf-stats.py', description=__doc__)
parser.add_argument('panel', metavar='panel', help='biallelic panel variants.')
parser.add_argument('pangenie', metavar='pangenie', help='PanGenie biallelic genotyped variants for all unrelated samples.')
parser.add_argument('pangenie_all', metavar='pangenie_all', help='PanGenie biallelic genotyped variants for all samples (including children).')
args = parser.parse_args()

# compute statistics for input panel VCF
panel_stats = {}

for line in gzip.open(args.panel, 'rt'):
	if line.startswith('#'):
		continue
	fields = line.strip().split()
	# require bi-allelic vcf with IDs
	assert len([a for a in fields[4].split(',')]) == 1
	info_fields = {k.split('=')[0] : k.split('=')[1] for k in fields[7].split(';') if '=' in k}
	assert 'ID' in info_fields
	var_id = info_fields['ID']
	allele_stats = compute_allele_statistics(fields)
	panel_stats[var_id] = allele_stats

sys.stderr.write('Done with panel.\n')


# compute statistics for all unrelated samples
pangenie_stats = {}
quals = [0,200]

for line in gzip.open(args.pangenie, 'rt'):
	if line.startswith('#'):
		continue
	fields = line.strip().split()
	# require bi-allelic vcf with IDs
	assert len([a for a in fields[4].split(',')]) == 1
	info_fields = {k.split('=')[0] : k.split('=')[1] for k in fields[7].split(';') if '=' in k}
	assert 'ID' in info_fields
	var_id = info_fields['ID']
	allele_stats = compute_allele_statistics(fields)
	genotype_stats, counts = compute_genotype_statistics(fields, quals)
	assert 'UK' in info_fields
	uk = info_fields['UK']
	pangenie_stats[var_id] = [allele_stats, genotype_stats, uk, counts]

sys.stderr.write('Done with unrelated.\n')

# compute statistics for all samples (unrelated + related)
pangenie_all_stats = {}

for line in gzip.open(args.pangenie_all, 'rt'):
	if line.startswith('#'):
		continue
	fields = line.strip().split()
	# require bi-allelic vcf with IDs
	assert len([a for a in fields[4].split(',')]) == 1
	info_fields = {k.split('=')[0] : k.split('=')[1] for k in fields[7].split(';') if '=' in k}
	assert 'ID' in info_fields
	var_id = info_fields['ID']
	allele_stats = compute_allele_statistics(fields)
	genotype_stats, counts = compute_genotype_statistics(fields, quals)
	assert 'UK' in info_fields
	uk = info_fields['UK']
	pangenie_all_stats[var_id] = [allele_stats, genotype_stats, uk, counts]

sys.stderr.write('Done with all.\n')

# print stats for all IDs in genotypes VCF
header = [ 	'variant_id',
		'panel_allele_freq',
		'panel_alternative_alleles',
		'panel_total_alleles',
		'panel_unknown_alleles',

		'pangenie-unrelated_allele_freq',
		'pangenie-unrelated_alternative_alleles',
		'pangenie-unrelated_total_alleles',
		'pangenie-unrelated_unknown_alleles',
		'pangenie-unrelated_heterozygosity',
		'pangenie-unrelated_heterozygous_genotypes',
		'pangenie-unrelated_total_genotypes',
		'pangenie-unrelated_unique_kmers',

		'pangenie-all_allele_freq',
		'pangenie-all_alternative_alleles',
		'pangenie-all_total_alleles',
		'pangenie-all_unknown_alleles',
		'pangenie-all_heterozygosity',
		'pangenie-all_heterozygous_genotypes',
		'pangenie-all_total_genotypes',
		'pangenie-all_unique_kmers'
	]

for q in quals:
	header.append('pangenie-unrelated_GQ>=' + str(q))
	header.append('pangenie-all_GQ>=' + str(q))

print('\t'.join(header))

assert len(pangenie_stats) == len(pangenie_all_stats)

for var_id in pangenie_stats:
	if not var_id in panel_stats:
		continue
		
	line = [	var_id,
			panel_stats[var_id].af,
			panel_stats[var_id].ac,
			panel_stats[var_id].an,
			panel_stats[var_id].untyped,

			pangenie_stats[var_id][0].af,
			pangenie_stats[var_id][0].ac,
			pangenie_stats[var_id][0].an,
			pangenie_stats[var_id][0].untyped,
			pangenie_stats[var_id][1].heterozygosity,
			pangenie_stats[var_id][1].het,
			pangenie_stats[var_id][1].total,
			pangenie_stats[var_id][2],

			pangenie_all_stats[var_id][0].af,
			pangenie_all_stats[var_id][0].ac,
			pangenie_all_stats[var_id][0].an,
			pangenie_all_stats[var_id][0].untyped,
			pangenie_all_stats[var_id][1].heterozygosity,
			pangenie_all_stats[var_id][1].het,
			pangenie_all_stats[var_id][1].total,
			pangenie_all_stats[var_id][2]

		]
		
	# add counts for GQs
	for q in quals:
		line.append(str(pangenie_stats[var_id][3][q]))
		line.append(str(pangenie_all_stats[var_id][3][q]))
	assert len(line) == len(header)
	print('\t'.join(line))
