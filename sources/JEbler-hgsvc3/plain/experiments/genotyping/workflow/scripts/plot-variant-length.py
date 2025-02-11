import matplotlib
matplotlib.use('pdf')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
import argparse
import sys
import gzip
from collections import defaultdict

from variantclassifier import determine_class_from_alleles


def extract_variant_from_line(line):
	splitted = line.split()
	alleles = [splitted[3]] + [s for s in splitted[4].split(',')]
	assert len(alleles) == 2
	is_symbolic = alleles[1][0] == '<'
	varlen = None
	vartype = None
	af = None
	info = {i.split('=')[0] : i.split('=')[1] for i in splitted[7].split(';') if "=" in i}
	if is_symbolic:
		vartype = alleles[1][1:-1]
		assert 'SVLEN' in info
		varlen = abs(int(info['SVLEN']))
	else:
		vartype = determine_class_from_alleles(alleles[0], alleles[1])
		varlen = max([len(a) for a in alleles])
		if (vartype == 'INS') or (vartype == 'DEL'):
			# insertions and deletions always have one base more
			varlen -= 1
		assert varlen >= 1
	assert  'AF' in info
	if info['AF'] == '.':
		info['AF'] = "0.0"
	return vartype, varlen, float(info['AF'])
	

def plot(numbers, labels, outname, af):
	with PdfPages(outname) as pdf:
		print(-1*np.logspace(1.6,8,200)[::-1])
		print(np.logspace(1.6,8,200))
		bins = np.append(-1*np.logspace(1.6,8,300)[::-1],  np.logspace(1.6,8,300))
#		bins = np.logspace(-8,-1.6, num=200) + np.logspace(1.6, 8, num=200)
		print(bins)
		plt.figure(figsize=(15,5))
		for n, label in zip(numbers, labels):
			y,binEdges=np.histogram(n[1]+n[2],bins=bins)
			bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
			plt.plot(bincenters,y,'-', label=label)
		plt.title('Common variants (AF > ' + str(af) + ')')
		plt.xscale('symlog')
		plt.yscale('symlog')
		plt.xlabel('allele length')
		plt.ylabel('count')
		plt.legend(loc='upper right')
		pdf.savefig()
		plt.close()



parser = argparse.ArgumentParser(prog='plot-variant-length.py')
parser.add_argument('-c', '--callsets', nargs='+', required=True, help="List of callsets to consider.")
parser.add_argument('-n', '--names', nargs='+', required=True, help="List of callset names.")
parser.add_argument('-o', '--outfile', required=True, help='name of the output file.')
parser.add_argument('-a', '--allele-freq', type=float, required=True, help='allele frequency threshold')
args = parser.parse_args()


callset_to_stats = defaultdict(lambda: [[],[],[]])
for callsetfile, name in zip(args.callsets, args.names):
	total = 0
	insertion = 0
	deletion = 0
	for line in gzip.open(callsetfile, 'rt'):
		if line.startswith('#'):
			continue
		vartype, varlen, af = extract_variant_from_line(line)
		if af <= args.allele_freq:
			continue
		total += 1
		callset_to_stats[name][0].append(varlen)
		if vartype == "DEL":
			deletion += 1
			callset_to_stats[name][1].append(-varlen)
		elif vartype in ['INS', 'DUP']:
			insertion += 1
			callset_to_stats[name][2].append(varlen)
	print(name)
	print('total:', total)
	print('insertion:', insertion)
	print('deletion:', deletion)
	

numbers = []
labels = []
for l,n in callset_to_stats.items():
	numbers.append(n)
	labels.append(l)

plot(numbers, labels, args.outfile, args.allele_freq)

