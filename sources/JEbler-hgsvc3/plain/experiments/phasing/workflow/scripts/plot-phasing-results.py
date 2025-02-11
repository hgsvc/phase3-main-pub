import argparse
import sys
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from collections import defaultdict
import statistics



def plot_switch_errors(switch_error_rates, sample_to_status):

	plot_samples_prev = None

	plt.figure(figsize=(9,6))

	for callset in switch_error_rates:

		unrelated_samples = []
		unrelated_values = []

		parent_samples = []
		parent_values = []

		child_samples = []
		child_values = []

		for sample in sorted(switch_error_rates[callset].keys()):
			value = switch_error_rates[callset][sample]
			if sample_to_status[sample] == "child":
				child_samples.append(sample)
				child_values.append(value)
			elif sample_to_status[sample] == "parent":
				parent_samples.append(sample)
				parent_values.append(value)
			else:
				unrelated_samples.append(sample)
				unrelated_values.append(value)


		plot_samples = [u + " (unrelated)" for u in unrelated_samples] + [p + " (parent)" for p in parent_samples] + [c + " (child)" for c in child_samples]
		plot_values = unrelated_values + parent_values + child_values


		# print stats
		print(callset)
		print("unrelated samples:")
		print("Min: " + str(min(unrelated_values)))
		print("Max: " + str(max(unrelated_values)))
		print("Mean: " + str(statistics.mean(unrelated_values)))
		print("Median: " + str(statistics.median(unrelated_values)))


		print("parent samples:")
		print("Min: " + str(min(parent_values)))
		print("Max: " + str(max(parent_values)))
		print("Mean: " + str(statistics.mean(parent_values)))
		print("Median: " + str(statistics.median(parent_values)))


		print("child samples:")
		print("Min: " + str(min(child_values)))
		print("Max: " + str(max(child_values)))
		print("Mean: " + str(statistics.mean(child_values)))
		print("Median: " + str(statistics.median(child_values)))



		if plot_samples_prev:
			assert plot_samples_prev == plot_samples
		plot_samples_prev = plot_samples


		x_values = [x for x in range(1, len(plot_samples) + 1)]
		plt.plot(x_values, plot_values, marker='o', label=callset)
		plt.xticks(x_values, plot_samples, rotation=90)
		
	plt.ylabel('Switch Error Rate [%]')
	plt.title('Comparison against ' + args.truthsetname + ' phasing')
	plt.legend()
	plt.tight_layout()
	plt.savefig(args.outname)


parser = argparse.ArgumentParser(prog='plot-phasing-results.py', description=__doc__)
parser.add_argument('-tsvfiles', metavar='TSVFILES', required=True, nargs='+', help='List of TSV files from whatshap compare.')
parser.add_argument('-truthsetname', metavar='TRUTHSETNAME', required=True, help='Name of the truth set for figure title.' )
parser.add_argument('-trios', metavar='TRIOS', required=True, help='PED file with trio information.')
parser.add_argument('-outname', metavar='OUTNAME', required=True, help='Name of the output PDF.')
args = parser.parse_args()


switch_error_rates = defaultdict(dict)

for tsv_file in args.tsvfiles:
	df = pd.read_csv(tsv_file, sep='\t')
	sample = df['#sample'][0]
	switch_errors = sum(df['all_switches'])
	variants = sum(df['all_assessed_pairs'])
	region = tsv_file.split('_')[-2]
	switch_error_rates[region][sample] = (switch_errors / variants) * 100.0


sample_to_status = defaultdict(lambda: "unrelated")

for trio in open(args.trios, 'r'):
	fields = trio.strip().split()

	# determine status of first sample
	if (fields[1] != "NA") or (fields[2] != "NA"):
		sample_to_status[fields[0]] = "child"

	# determine status of second and third sample (if present)
	for sample in [fields[1], fields[2]]:
		if sample == "NA":
			continue
		if sample_to_status[sample] != "child":
			sample_to_status[sample] = "parent"

plot_switch_errors(switch_error_rates, sample_to_status)
