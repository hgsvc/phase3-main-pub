import sys
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import argparse


parser = argparse.ArgumentParser(prog='plot-qv.py', description=__doc__)
parser.add_argument('outname', metavar='OUTNAME', help='Name of the output PDF.')
parser.add_argument('-assembly', metavar='ASSEMBLYQVs', required=False, default='')
args = parser.parse_args()

merged = {}
single = {}
genome = {}
assembly = {}

samples = set([])

for filename in sys.stdin:
	fields = filename.strip().split('/')[-1].split('_')
	mode = fields[0]
	sample = fields[1]
	samples.add(sample)
	haplotype = fields[2].split('.')[0]

	for line in open(filename.strip(), 'r'):
		if line.startswith('QV'):
			adjusted_qv = line.strip().split()[-1]
			assert mode in ['merged', 'single', 'genome']
			if mode == 'merged':
				merged[sample + '_' + haplotype] = float(adjusted_qv)
			elif mode == 'genome':
				genome[sample + '_' + haplotype] = float(adjusted_qv)
			else:
				single[sample + '_' + haplotype] = float(adjusted_qv)


if args.assembly:
	for line in open(args.assembly, 'r'):
		fields = line.strip().split()
		if line.startswith('hap'):
			continue
		if fields[0] in ['H0', 'wg']:
			continue
		sample = fields[-2].split('.')[0]
		assert fields[0] in ['H1', 'H2']
		haplotype = 'hap1' if fields[0] == 'H1' else 'hap2'
		assembly[sample + '_' + haplotype] = float(fields[-1])


# plot results
y_merged = []
y_single = []
y_genome = []
y_assembly = []

labels = []

for sample in sorted(list(samples)):
	for haplotype in ['hap1', 'hap2']:
		if sample == "HG00732":
			continue
		labels.append(sample + '_' + haplotype)
		y_merged.append(merged[sample + '_' + haplotype])
		y_single.append(single[sample + '_' + haplotype])
		y_genome.append(genome[sample + '_' + haplotype])
		if args.assembly:
			y_assembly.append(assembly[sample + '_' + haplotype])


with PdfPages(args.outname) as pdf:
	x = [i*2 for i in range(len(labels))]

	plt.figure(figsize=(20,10))
	plt.plot(x, y_merged, marker = 's', label = 'merged calls', color = '#377eb8')
	plt.plot(x, y_single, marker = 'd', label = 'single sample calls', color = '#ff7f00')
	plt.plot(x, y_genome, marker= '^', label = 'reference sequence', color = '#4daf4a')
	if args.assembly:
		plt.plot(x, y_assembly, marker= 'o', label = 'assembly', color = '#f781bf')
	plt.xticks(x, labels, rotation = 'vertical')
	plt.legend()
	plt.ylabel('QV')
	plt.tight_layout()
	pdf.savefig()
	plt.close()

	plt.figure()
	if args.assembly:
		plt.boxplot([y_assembly, y_single, y_merged, y_genome])
		plt.xticks([1,2,3,4], ['assembly', 'single sample calls', 'merged calls', 'reference sequence'], rotation='vertical')
	else:
		plt.boxplot([y_single, y_merged, y_genome])
		plt.xticks([1,2,3], ['single sample calls', 'merged_calls', 'reference sequence'], rotation='vertical')
	plt.ylabel('QV')
	plt.tight_layout()
	pdf.savefig()
