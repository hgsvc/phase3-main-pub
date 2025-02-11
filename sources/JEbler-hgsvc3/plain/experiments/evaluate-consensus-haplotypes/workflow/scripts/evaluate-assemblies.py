import sys
import argparse
from collections import defaultdict
import gzip
import math

def compute_length(filename):
	total = 0
	for line in gzip.open(filename, 'rt'):
		if line.startswith('>'):
			continue
		total += len(line.strip())
	return total



def analyze_gt_cons(allele_h, allele_assemb1, allele_assemb2):
	observed = [allele_assemb1, allele_assemb2, allele_h]
	if observed == ['1','0','0']:
		return 'FN_UNKNOWN'
	elif observed == ['0','1','0']:
		return 'FN_UNKNOWN'
	elif observed == ['1','1','0']:
		return 'FN_BOTH'
	elif observed == ['0','0','1']:
		return 'FP_H2'
	elif observed == ['1','0','1']:
		return 'TP'
	elif observed == ['0','1','1']:
		return 'TP'
	elif observed == ['1','1','1']:
		return 'FN_H1'
	
	elif observed == ['1', '1', '.']:
		return 'FN_UNKNOWN'
	elif observed == ['1', '.', '0']:
		return 'FN_UNKNOWN'
	elif observed == ['.', '1', '0']:
		return 'FN_UNKNOWN'
	else:
		assert observed in [ ['1', '.', '1'], ['1', '.', '.'], ['.', '1', '.'], ['.', '.', '1'], ['.', '0', '1'], ['1', '0', '.'], ['.', '1', '1'], ['0', '.', '1'], ['0', '1', '.']]
		return 'UNKNOWN'



def analyze_gt_assemb(allele_assemb, allele_h1, allele_h2):
	observed = [allele_assemb, allele_h1, allele_h2]
	if observed == ['1','0','0']:
		return 'FN_UNKNOWN'
	elif observed == ['0','1','0']:
		return 'FP_H1'
	elif observed == ['1','1','0']:
		return 'TP'
	elif observed == ['0','0','1']:
		return 'FP_H2'
	elif observed == ['1','0','1']:
		return 'TP'
	elif observed == ['0','1','1']:
		return 'FP_BOTH'
	elif observed == ['1','1','1']:
		return 'FP_UNKNOWN'
	elif observed == ['.', '1', '1']:
		return 'FP_UNKNOWN'
	elif observed == ['0', '.', '1']:
		return 'FP_H2'
	elif observed == ['0', '1', '.']:
		return 'FP_H1'
	else:
		assert observed in [['1', '1', '.'], ['1', '.', '1'], ['1', '.', '.'], ['.', '1', '.'], ['.', '.', '1'], ['.', '0', '1'], ['1', '0', '.'], ['1', '.', '0'], ['.', '1', '0']]
		return 'UNKNOWN'





if __name__ == '__main__':

	parser = argparse.ArgumentParser(prog='evaluate-assemblies.py', description=__doc__)
	parser.add_argument('--name', required = True, metavar='NAME', help='Name of the reference haplotype used (to print).')
	parser.add_argument('--qv-h1', required = False, metavar='QV1', default="nan", help="QV value computed for first haplotype.")
	parser.add_argument('--qv-h2', required = False, metavar='QV2', default="nan", help="QV value computed for second haplotype.")
	parser.add_argument('--hap1', required = True, metavar='HAP1', help="FASTA first haplotype.")
	parser.add_argument('--hap2', required = True, metavar='HAP2', help="FASTA second haplotype.")
	parser.add_argument('--errors', required = True, metavar='ERROR', help="Write errors into this BED file.")
	parser.add_argument('--all', required = True, metavar='ALL', help="Write all variants into this BED file.")
	parser.add_argument('--reference', required = True, choices = ['assembly', 'consensus'], help="Which sequence was used as a reference.")
	args = parser.parse_args()


	bp_changes_error = 0
	stats = {}
	stats['FN_UNKNOWN'] = 0
	stats['FP_UNKNOWN'] = 0
	stats['FP_H1'] = 0
	stats['FN_H1'] = 0
	stats['TP'] = 0
	stats['FP_H2'] = 0
	stats['FN_H2'] = 0
	stats['FP_BOTH'] = 0
	stats['FN_BOTH'] = 0
	stats['FP_UNKNOWN'] = 0
	stats['FN_UNKNOWN'] = 0
	stats['UNKNOWN'] = 0
	total = 0
	skipped = 0

	len_h1 = compute_length(args.hap1)
	len_h2 = compute_length(args.hap2)

	with open(args.errors, 'w') as outfile_err, open(args.all, 'w') as outfile_all:

		for line in sys.stdin:
			if line.startswith('#'):
				continue
			fields = line.strip().split()
			if fields[6] != 'PASS':
				continue
			if 'chrY' in line:
				continue
			ref_allele = fields[3]
			alt_alleles = fields[4].split(',')
			assert len(alt_alleles) == 1
			alt_allele = alt_alleles[0]
			info_fields = {f.split('=')[0] : f.split('=')[1] for f in fields[7].split(';') if '=' in f}
			assert 'ID' in info_fields
			var_id = info_fields['ID']
			bp_change = 1 if 'SNV' in var_id else abs(int(info_fields['SVLEN']))
			gt = fields[9].strip().split('|')
			label = analyze_gt_assemb(gt[2], gt[0], gt[1]) if (args.reference == "assembly") else analyze_gt_cons(gt[2], gt[0], gt[1])
			end = int(fields[1]) + len(fields[3])
			if label:
				outfile_all.write('\t'.join([fields[0], fields[1], str(end), str(bp_change), info_fields['ID'], label]) + '\n')
				stats[label] += 1
				if label not in ['TP', 'UNKNOWN']:
					bp_changes_error += bp_change
					outfile_err.write('\t'.join([fields[0], fields[1], str(end), str(bp_change), info_fields['ID'], label]) + '\n')
					if label in ['FP_BOTH', 'FN_BOTH']:
						bp_changes_error += bp_change
			else:
				skipped += 1
			total += 1

		# print statistics
		keys = sorted([k for k in stats.keys()])
		score = bp_changes_error / (len_h1 + len_h2)
		print('\t'.join(['ref_HT', 'total_vars', 'skipped_vars'] + keys + ['total_length', 'bp_changes', 'score', 'score_phred', 'QV_h1', 'QV_h2', 'QV_avg']))
		avg_qv = 'nan' if (args.qv_h1 == 'nan') or (args.qv_h2 == 'nan') else (float(args.qv_h1) + float(args.qv_h2)) / 2
		stats_to_print = [args.name, str(total), str(skipped)] + [str(stats[k]) for k in keys] + [str(len_h1 + len_h2), str(bp_changes_error), str(score), str(-10*math.log10(score)), args.qv_h1, args.qv_h2, str(avg_qv)]
		print('\t'.join(stats_to_print))
