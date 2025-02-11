import sys
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import statistics

phred_scores = []
nr_svs = []

for line in sys.stdin:
	if line.startswith('ref_HT'):
		continue
	fields = line.strip().split()
	if fields[16] != 'inf':
		phred_scores.append(float(fields[16]))
		nr_svs.append(int(fields[18]))


with PdfPages(sys.argv[1]) as pdf:

	# plot phred score historgram
	plt.figure()
	plt.title(sys.argv[2])
	plt.hist(phred_scores, bins = range(int(min(phred_scores)) - 1, int(max(phred_scores)) + 2, 1))
	plt.axvline(statistics.median(phred_scores), color='k', linestyle='dashed', linewidth=1)

	plt.xlim([0,65])
	plt.ylim([0,500])
	plt.xlabel('Phred scores')
	plt.ylabel('Count')
	pdf.savefig()

	print('Phred scores:')
	print('min: ' + str(min(phred_scores)))
	print('max: ' + str(max(phred_scores)))
	print('mean: ' + str(statistics.mean(phred_scores)))
	print('median: ' + str(statistics.median(phred_scores)))

	# plot histogram on the number of SVs by interval
	plt.figure()
	plt.title(sys.argv[2])
	plt.hist(nr_svs, bins = range(int(min(nr_svs)) - 1, int(max(nr_svs)) + 2, 1))
	plt.axvline(statistics.median(nr_svs), color='k', linestyle='dashed', linewidth=1)

	plt.xlim([0,150])
	plt.ylim([0,400])
	plt.xlabel('Number of variants (> 20bp)')
	plt.ylabel('Count')
	pdf.savefig()

	print('Number of variants per interval:')
	print('min: ' + str(min(nr_svs)))
	print('max: ' + str(max(nr_svs)))
	print('mean: ' + str(statistics.mean(nr_svs)))
	print('median: ' + str(statistics.median(nr_svs)))


