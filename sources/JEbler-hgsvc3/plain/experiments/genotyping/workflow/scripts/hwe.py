import sys
from collections import namedtuple
from collections import defaultdict
import argparse
import vcf
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from enum import Enum
import numpy as np


Position = namedtuple('Position', 'chromosome start')

def benjamini_hochberg(p_values, alpha):
	n_tests = len(p_values)
	# sort in ascending order
	sorted_pvals = sorted(p_values, key=lambda x: x[1])
	values = [ ((i+1)/n_tests)*alpha for i in range(n_tests)]
	# determine largest index with p_val < val
	index = -1
	for i in range(n_tests):
		if sorted_pvals[i][1] <= values[i]:
			index = i
	# return significant p_values
	return sorted_pvals[0:index+1]


parser = argparse.ArgumentParser(prog='hwe.py', description=__doc__)
parser.add_argument('hwe', metavar='HWE', nargs='+', help="hwe files from VCFtools for each variant category.")
parser.add_argument('--labels', metavar='LABELS', required=True, nargs='+', help='labels in order of first argument.')
parser.add_argument('--outname', metavar='OUTNAME', help='prefix of output file names.', default='output')
args = parser.parse_args()

filenames = args.hwe
labels = args.labels

assert len(filenames) == len(labels)

p_values = []
type_to_count = defaultdict(lambda: 0)
for filename,label in zip(filenames, labels):
	nr_lines = 0
	for line in open(filename, 'r'):
		if line.startswith('CHR'):
			continue
		nr_lines += 1
		fields = line.split()
		genotype_counts = [int(i) for i in fields[2].split('/')]
		n_samples = sum(genotype_counts)
		het = (genotype_counts[1] / n_samples)
		af = (genotype_counts[1] + 2*genotype_counts[2]) / (2*n_samples)

		if af > 0.0:
			p_value = float(fields[5])
			position = Position(fields[0], int(fields[1]))
			vartype = label
			type_to_count[vartype] += 1
			p_values.append((position, p_value, vartype))
	sys.stderr.write("Read " + str(nr_lines) + " " + filename + '. \n')


# compute sigificant p-values
significant = benjamini_hochberg(p_values, 0.05)
# determine significants per type
significant_per_type = defaultdict(lambda: 0)
for p in significant:
	significant_per_type[p[2]] += 1

# print statistics
with open(args.outname + '_hwe.tsv', 'w') as outhwe:
	# header
	header = '\t'.join([
		'type',
		'total',
		'not_significant',
		'not_significant[%]']) + '\n'

	outhwe.write(header)

	total = 0
	total_hwe = 0

	for label in labels:
		nr_total = type_to_count[label]
		nr_hwe = nr_total - significant_per_type[label]
		percentage_hwe = (float(nr_hwe) / max(1.0, float(nr_total))) * 100.0
		total += nr_total
		total_hwe += nr_hwe
		line = [label,
			str(nr_total),
			str(nr_hwe),
			str(percentage_hwe)]
		outhwe.write('\t'.join(line) + '\n')

	# overall statistics
	assert total == len(p_values)
	percent_all = (float(total_hwe) / float(total)) * 100.0
	line = [
		'all',
		str(total),
		str(total_hwe),
		str(percent_all)
	]
	outhwe.write('\t'.join(line) + '\n')
