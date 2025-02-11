import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import argparse


def plot_concordances_all(files, outname, sources):
	var_to_name = {
		'snp' : 'SNPs',
		'indels': 'indels (1-49bp)',
		'large-insertion': 'SV insertions (>=50bp)',
		'large-deletion': 'SV deletions (>=50bp)',
		'large-complex': 'SV complex (>=50bp)'
		}

	colors = ['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3', '#808080', '#f7ce37', '#bc202c', '#251188']
	type_to_file = {}
	n_rows = 1
	n_cols = 5
	plt.figure(figsize=(65,5))
	for source in sources:
		for f in files:
			if not ('_' + source + '_') in f:
				continue
			vartype = f.split('_')[-1][:-4]
			type_to_file[(source, vartype)] = f
	variants = ['snp', 'indels', 'large-deletion', 'large-insertion', 'large-complex']
	plot_index = 1
	for var in variants:
		plt.subplot(n_rows, n_cols, plot_index)
		plot_index += 1
		all_samples = []
		x_values = []
		is_first = True
		for i,source in enumerate(sources):
			print('source', source)
			samples = []
			concordances = []
			if not (source, var) in type_to_file:
				continue
			for line in open(type_to_file[(source,var)], 'r'):
				if line.startswith('sample'):
					continue
				fields = line.split()
				samples.append(fields[0])
				concordances.append(float(fields[1]))
			if is_first:
				all_samples = samples
			else:
				assert all_samples == samples
			is_first = False
			x_values = [i*6 for i in range(len(samples))]
			plt.plot(x_values, concordances, label=source, color=colors[i], marker='o')
		plt.ylim(0.0,100)
		plt.title(var_to_name[var])
		plt.xticks(x_values, all_samples, rotation='vertical')
		plt.ylabel('weighted genotype concordance [%]')
		plt.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.3)
		plt.tight_layout()
	# create legend
	handles = []
	labels = []
	for i, source in enumerate(sources):
		label = source
		line = matplotlib.lines.Line2D([],[], color=colors[i], markersize=100, linewidth=2.0, linestyle='-', label=label)
		handles.append(line)
		labels.append(label)
	plt.figlegend(handles, labels)
	plt.savefig(outname)



def plot_untyped_all(files, outname, sources):
	var_to_name = {
		'snp' : 'SNPs',
		'indels': 'indels (1-49bp)',
		'large-insertion': 'SV insertions (>=50bp)',
		'large-deletion': 'SV deletions (>=50bp)',
		'large-complex': 'SV complex (>=50bp)'
		}

	colors = ['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3', '#808080', '#f7ce37', '#bc202c', '#251188']
	type_to_file = {}
	n_rows = 1
	n_cols = 5
	plt.figure(figsize=(65,5))
	for source in sources:
		for f in files:
			if not ('_' + source + '_') in f:
				continue
			vartype = f.split('_')[-1][:-4]
			type_to_file[(source, vartype)] = f
	variants = ['snp', 'indels', 'large-deletion', 'large-insertion', 'large-complex']
	plot_index = 1
	for var in variants:
		plt.subplot(n_rows, n_cols, plot_index)
		plot_index += 1
		all_samples = []
		x_values = []
		is_first = True
		for i,source in enumerate(sources):
			print('source', source)
			samples = []
			untyped = []
			if not (source, var) in type_to_file:
				continue
			for line in open(type_to_file[(source,var)], 'r'):
				if line.startswith('sample'):
					continue
				fields = line.split()
				samples.append(fields[0])
				untyped.append(100.0 - float(fields[2]))
			if is_first:
				all_samples = samples
			else:
				assert all_samples == samples
			is_first = False
			x_values = [i*6 for i in range(len(samples))]
			plt.plot(x_values, untyped, label=source, color=colors[i], marker='o')
		plt.ylim(0.0, 7.0)
		plt.title(var_to_name[var])
		plt.xticks(x_values, all_samples, rotation='vertical')
		plt.ylabel('untyped variants [%]')
		plt.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.3)
		plt.tight_layout()
	# create legend
	handles = []
	labels = []
	for i, source in enumerate(sources):
		label = source
		line = matplotlib.lines.Line2D([],[], color=colors[i], markersize=100, linewidth=2.0, linestyle='-', label=label)
		handles.append(line)
		labels.append(label)
	plt.figlegend(handles, labels)
	plt.savefig(outname)



def plot_fscores_all(files, outname, sources):	
	var_to_name = {
		'snp' : 'SNPs',
		'indels': 'indels (1-49bp)',
		'large-insertion': 'SV insertions (>=50bp)',
		'large-deletion': 'SV deletions (>=50bp)',
		'large-complex': 'SV complex (>=50bp)'
		}

	colors = ['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3', '#808080', '#f7ce37', '#bc202c', '#251188']
	type_to_file = {}
	n_rows = 1
	n_cols = 5
	plt.figure(figsize=(65,5))
	for source in sources:
		for f in files:
			if not ('_' + source + '_') in f:
				continue
			vartype = f.split('_')[-1][:-4]
			type_to_file[(source, vartype)] = f
	variants = ['snp', 'indels', 'large-deletion', 'large-insertion', 'large-complex']
	plot_index = 1
	for var in variants:
		plt.subplot(n_rows, n_cols, plot_index)
		plot_index += 1
		all_samples = []
		is_first = True
		for i,source in enumerate(sources):
			print('source', source)
			samples = []
			fscores = []
			if not (source,var) in type_to_file:
				continue
			for line in open(type_to_file[(source,var)], 'r'):
				if line.startswith('sample'):
					continue
				fields = line.split()
				samples.append(fields[0])
				fscores.append(float(fields[3]))
			if is_first:
				all_samples = samples
			else:
				assert all_samples == samples
			is_first = False
			x_values = [i*6 for i in range(len(samples))]
			plt.plot(x_values, fscores, label=source, color=colors[i], marker='o')
		plt.title(var_to_name[var])
		plt.xticks(x_values, all_samples, rotation='vertical')
		plt.ylabel('adjusted F-score [%]')
		plt.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.3)
		plt.tight_layout()
	# create legend
	handles = []
	labels = []
	for i, source in enumerate(sources):
		label = source
		line = matplotlib.lines.Line2D([],[], color=colors[i], markersize=100, linewidth=2.0, linestyle='-', label=label)
		handles.append(line)
		labels.append(label)
	plt.figlegend(handles, labels)
	plt.savefig(outname)



def plot_concordance_vs_untyped(files, outname, sources):
	var_to_name = {
		'snp' : 'SNPs',
		'indels': 'indels (1-49bp)',
		'large-insertion': 'SV insertions (>=50bp)',
		'large-deletion': 'SV deletions (>=50bp)',
		'large-complex': 'SV complex (>=50bp)'
		}
	colors = ['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3', '#808080', '#f7ce37', '#bc202c', '#251188']
	type_to_file = {}
	n_rows = 5 # one row per variant type
	plt.figure(figsize=(80,30))
	for source in sources:
		for f in files:
			if not ('_' + source + '_') in f:
				continue
			vartype = f.split('_')[-1][:-4]
			type_to_file[(source, vartype)] = f
	variants = ['snp', 'indels', 'large-deletion', 'large-insertion', 'large-complex']
	data = {}
	all_samples = []
	is_first = True
	plot_index = 1
	for var in variants:
		for i, source in enumerate(sources):
			if not (source, var) in type_to_file:
				continue	
			for line in open(type_to_file[(source,var)], 'r'):
				if line.startswith('sample'):
					continue
				fields = line.split()
				sample = fields[0]
				if is_first:
					all_samples.append(sample)
				concordance = float(fields[1])
				typed = float(fields[2])
				data[(var, source, sample)] = (concordance, typed)
			is_first = False
	for var in variants:
		for sample in all_samples:
			plt.subplot(n_rows, len(all_samples), plot_index)
			for i,source in enumerate(sources):
				if (var, source, sample) in data:
					plt.plot(data[(var, source, sample)][1], data[(var, source, sample)][0],  marker='o', label=source, color=colors[i])
			plt.title(var_to_name[var] + ', ' + sample)
			plt.ylabel('weighted genotype concordance [%]')
			plt.xlabel('typed variants [%]')
			plt.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.3)
			plt.tight_layout()
			plot_index += 1
	# create legend
	handles = []
	labels = []
	for i, source in enumerate(sources):
		label = source
		line = matplotlib.lines.Line2D([],[], markersize=100, linewidth=2.0, linestyle='-', label=label, color=colors[i])
		handles.append(line)
		labels.append(label)
	plt.figlegend(handles, labels)
	plt.savefig(outname)



parser = argparse.ArgumentParser(prog='plot-results.py', description="Plot concordances or precision/recall statistics.")
parser.add_argument('-files', metavar='FILES', nargs='+', help='files with results per sample.')
parser.add_argument('-outname', metavar='OUTNAME', required=True, help='Name of the output file.')
parser.add_argument('-sources',  metavar='SOURCES', nargs='+', help='regions concordance_all')
parser.add_argument('-metric', metavar='METRIC', required=True, choices=['concordance', 'precision-recall-typable', 'untyped', 'concordance-vs-untyped'], help='Which metric to use for plotting.')
args = parser.parse_args()

if args.metric == 'concordance':
	plot_concordances_all(args.files, args.outname, args.sources)
elif args.metric == 'untyped':
	plot_untyped_all(args.files, args.outname, args.sources)
elif args.metric == 'precision-recall-typable':
	plot_fscores_all(args.files, args.outname, args.sources)
else:
	plot_concordance_vs_untyped(args.files, args.outname, args.sources)
