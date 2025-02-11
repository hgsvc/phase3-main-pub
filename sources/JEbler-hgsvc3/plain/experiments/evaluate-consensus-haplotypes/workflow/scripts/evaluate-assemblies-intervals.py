import sys
import argparse
from collections import defaultdict
import gzip
import math
import statistics

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


def parse_chromosome_length(line):
	assert "ID=" in line
	assert "length=" in line
	chrom = line.strip().split(',')[0].split('=')[-1]
	length = int(line.strip().split(',')[1].split('=')[-1])
	return chrom, length


class Statistics:
	def __init__(self, start, end, length):
		self.bp_changes_error = 0
		self.nr_svs = 0
		self.stats = {}
		self.stats['FN_UNKNOWN'] = 0
		self.stats['FP_UNKNOWN'] = 0
		self.stats['FP_H1'] = 0
		self.stats['FN_H1'] = 0
		self.stats['TP'] = 0
		self.stats['FP_H2'] = 0
		self.stats['FN_H2'] = 0
		self.stats['FP_BOTH'] = 0
		self.stats['FN_BOTH'] = 0
		self.stats['FP_UNKNOWN'] = 0
		self.stats['FN_UNKNOWN'] = 0
		self.stats['UNKNOWN'] = 0
		self.total = 0
		self.skipped = 0

		self.start = start
		self.end = end
		self.length = length

	def compute_score(self):
		return self.bp_changes_error / (self.length)

	def compute_phred_score(self):
		score = self.compute_score()
		return -10 * math.log10(score) if score > 0 else 'inf'

	def to_string(self, name, qv_h1, qv_h2, skip_header, interval_score = 'nan'):
		keys = sorted([k for k in self.stats.keys()])
		score = self.compute_score()
		stats_to_print = ""
		if not skip_header:
			stats_to_print = '\t'.join(['ref_HT', 'total_vars', 'skipped_vars'] + keys + ['total_length', 'bp_changes', 'score', 'score_phred', 'score_phred_interval', 'nr_sv_err', 'QV_h1', 'QV_h2', 'QV_avg']) + '\n'
		avg_qv = 'nan' if (qv_h1 == 'nan') or (qv_h2 == 'nan') else (float(qv_h1) + float(qv_h2)) / 2
		stats_to_print += '\t'.join([name, str(self.total), str(self.skipped)] + [str(self.stats[k]) for k in keys] + [str(self.length), str(self.bp_changes_error), str(score), str(self.compute_phred_score()), str(interval_score), str(self.nr_svs), qv_h1, qv_h2, str(avg_qv)])
		return stats_to_print



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
	parser.add_argument('--intervals', required = True, metavar='INTERVAL', help="Write interval based scores into this file.")
	args = parser.parse_args()

	chrom_to_length = {}
	len_h1 = compute_length(args.hap1)
	len_h2 = compute_length(args.hap2)
	interval_length = 1000000
	
	total_stats = Statistics(0, 0, len_h1+len_h2)
	interval_stats = defaultdict(list)


	with open(args.errors, 'w') as outfile_err, open(args.all, 'w') as outfile_all, open(args.intervals, 'w') as outfile_intervals:

		for line in sys.stdin:
			if line.startswith('##'):
				if "contig=<ID=" in line:
					chrom, length = parse_chromosome_length(line)
					chrom_to_length[chrom] = length
				continue
			if line.startswith('#'):
				# all header lines read. Initialize intervals
				for chrom, length in chrom_to_length.items():
					for interval_id in range(0, (length // interval_length) + 1):
						interval_start = interval_id*interval_length
						interval_end = min((interval_id + 1) * interval_length, length)
						interval_stats[chrom].append(Statistics(interval_start, interval_end, 2*(interval_end - interval_start)))
				continue
			fields = line.strip().split()
			if fields[6] != 'PASS':
				continue
			if 'chrY' in line:
				continue

			chrom = fields[0]
			ref_allele = fields[3]
			alt_alleles = fields[4].split(',')
			assert len(alt_alleles) == 1
			alt_allele = alt_alleles[0]
			info_fields = {f.split('=')[0] : f.split('=')[1] for f in fields[7].split(';') if '=' in f}
			assert 'ID' in info_fields
			var_id = info_fields['ID']
			bp_change = 1 if 'SNV' in var_id else abs(int(info_fields['SVLEN']))
			is_sv = False if 'SNV' in var_id else abs(int(info_fields['SVLEN'])) > 20
			gt = fields[9].strip().split('|')
			label = analyze_gt_assemb(gt[2], gt[0], gt[1]) if (args.reference == "assembly") else analyze_gt_cons(gt[2], gt[0], gt[1])
			end = int(fields[1]) + len(fields[3])
			interval_id = int(fields[1]) // interval_length
			if label:
				outfile_all.write('\t'.join([fields[0], fields[1], str(end), str(bp_change), info_fields['ID'], label]) + '\n')
				total_stats.stats[label] += 1
				interval_stats[chrom][interval_id].stats[label] += 1
				
				if label not in ['TP', 'UNKNOWN']:
					total_stats.bp_changes_error += bp_change
					total_stats.nr_svs += 1 if is_sv else 0
					interval_stats[chrom][interval_id].bp_changes_error += bp_change
					interval_stats[chrom][interval_id].nr_svs += 1 if is_sv else 0

					outfile_err.write('\t'.join([fields[0], fields[1], str(end), str(bp_change), info_fields['ID'], label]) + '\n')
					if label in ['FP_BOTH', 'FN_BOTH']:
						total_stats.bp_changes_error += bp_change
						total_stats.nr_svs += 1 if is_sv else 0
						interval_stats[chrom][interval_id].bp_changes_error += bp_change
						interval_stats[chrom][interval_id].nr_svs += 1 if is_sv else 0
			else:
				total_stats.skipped += 1
				interval_stats[chrom][interval_id].skipped += 1
			total_stats.total += 1
			interval_stats[chrom][interval_id].total += 1

		# print statistics
		skip_header = False
		interval_scores = []
		for chrom, stats in interval_stats.items():
			for s in stats:
				interval_name = args.name + '_' + chrom + ':' + str(s.start) + '-' + str(s.end)
				outfile_intervals.write(s.to_string(interval_name, 'nan', 'nan', skip_header) + '\n')
				skip_header = True
				phred = s.compute_phred_score()
				if phred != 'inf':
					interval_scores.append(phred)

		print(total_stats.to_string(args.name, args.qv_h1, args.qv_h2, False, statistics.median(interval_scores)))

