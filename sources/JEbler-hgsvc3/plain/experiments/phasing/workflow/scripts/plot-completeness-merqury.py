import sys
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import argparse
from collections import defaultdict
import statistics

parser = argparse.ArgumentParser(prog='plot-completeness-merqury.py', description=__doc__)
parser.add_argument('outname', metavar='OUTNAME', help='Name of the output PDF.')
args = parser.parse_args()


samples = set([])
comp_values = defaultdict(lambda: {})
callset_samples = set([])
assembly_samples = set([])

for filename in sys.stdin:
	fields = filename.strip().split('/')[-1].split('_')
	mode = filename.strip().split('/')[-5]
	callset = fields[0] + ( "" if mode == "assigned" else "-" + mode)
	callset += " (" + fields[1] + ")"
	sample = fields[2]
	haplotype = fields[3].split('.')[0]
	assert haplotype in ['hap1', 'hap2']
	samples.add(sample)
	callset_samples.add(sample)

	for line in open(filename.strip(), 'r'):
		completeness = line.strip().split()[-1]
		comp_values[callset][sample + '_' + haplotype] = float(completeness)


y_values = defaultdict(lambda: [])
labels = []

ordered_samples = sorted([s for s in samples if not s in assembly_samples]) + sorted([s for s in samples if s in assembly_samples])

for sample in ordered_samples:
	labels.append(sample + '_hap1')
	labels.append(sample + '_hap2')
	for callset in comp_values.keys():
		for haplotype in ['hap1', 'hap2']:
			key_str = sample + '_' + haplotype
			value = comp_values[callset][key_str] if key_str in comp_values[callset] else None
			y_values[callset].append(value)




	
with PdfPages(args.outname) as pdf:
	x = [i*2 for i in range(len(labels))]

	plt.figure(figsize=(30,10))
	for callset in y_values.keys():
		plt.plot(x, y_values[callset], marker = 'o', label = callset)
	plt.xticks(x, labels, rotation = 'vertical')
	plt.legend()
	plt.ylabel('completeness')
	plt.tight_layout()
	pdf.savefig()
	plt.close()


	values = [ [y for y in y_values[k] if y is not None] for k in sorted(y_values.keys()) ]
	x_labels = [k for k in sorted(y_values.keys())]

	plt.figure(figsize=(10,15))
	plt.boxplot(values)
	plt.xticks([i+1 for i in range(len(x_labels))], x_labels, rotation='vertical')
	plt.ylabel('completeness')
	plt.tight_layout()
	pdf.savefig()


	for v,k in zip(values, x_labels):
		print(k)
		print("Median: " +  str(statistics.median(v)))
		print("Min: " + str(min(v)))
		print("Max: " + str(max(v)))





