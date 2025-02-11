import sys
from collections import defaultdict
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

def is_fully_called(interval_length, var_length):
	if interval_length * var_length < 0:
		# signs are different
		return False

	if abs(interval_length)*0.5 < abs(var_length) < abs(interval_length)*1.5:
		return True
	else:
		return False

def is_partially_called(interval_length, var_length):
	if interval_length * var_length < 0:
		# signs are different
		return False

	if abs(var_length) > 0.2*abs(interval_length):
		return True
	else:
		return False	
	


def plot(interval_to_samples, interval_to_changes, interval_to_status, min_len, max_len, outprefix, missed_mc):

	plt.figure(figsize=(30.0,10.0))
	plt.rcParams.update({'font.size': 12})
	x_ticks = []
	interval_id = 0

	for interval in interval_to_samples:
	
		interval_len = int(interval.split('_')[-1])
		if (abs(interval_len) < min_len) or (abs(interval_len) > max_len):
			continue

		all_samples_to_plot = []
		carrier_samples_to_plot = []
		x_ticks.append(interval)

		for i in range(len(all_samples)):
			sample = all_samples[i]
			all_samples_to_plot.append(interval_to_changes[interval][i][0])
			all_samples_to_plot.append(interval_to_changes[interval][i][1])
			if sample in interval_to_samples[interval]:
				carrier_samples_to_plot.append(interval_to_changes[interval][i][0])
				carrier_samples_to_plot.append(interval_to_changes[interval][i][1])

		plt.scatter([interval_id*2] * len(all_samples_to_plot), all_samples_to_plot, color="grey", marker='o', alpha=0.3)
		plt.scatter([interval_id*2] * len(carrier_samples_to_plot), carrier_samples_to_plot, marker='o', color="blue", alpha=0.6)
		plt.scatter([interval_id*2], [ interval_len ], marker = "X", color="red")
		interval_id += 1


	tick_colors = []
	for t in x_ticks:
		if interval_to_status[t] == "FULLY_CALLED":
			tick_colors.append("green")
		elif interval_to_status[t] == "PARTIALLY_CALLED":
			tick_colors.append("orange")
		elif interval_to_status[t] == "NO_MATCHING_CALLS":
			tick_colors.append("black")
		else:
			tick_colors.append("darkgrey")

	marked_x_ticks = []
	for x in x_ticks:
		if missed_mc[x]:
			marked_x_ticks.append(x + '(*)')
		else:
			marked_x_ticks.append(x)
	x_ticks = marked_x_ticks

	# plot a line at zero
	line = plt.plot([i*2 for i in range(len(x_ticks))], [0.0] * len(x_ticks), color="grey")
	plt.setp(line, linewidth=0.4)   

	plt.xticks([i*2 for i in range(len(x_ticks))], x_ticks, rotation=90)
	
	for i,c in enumerate(tick_colors):
        	plt.gca().get_xticklabels()[i].set_color(c)
	plt.margins(x=0.01)
	plt.ylim(-max_len*2, max_len*2)
	plt.title("gene-interruptive SVs (" + str(min_len) + "-" + str(max_len) + "bp)")
	plt.ylabel("length change wrt. reference [bp]")
	plt.tight_layout()
	plt.savefig(outprefix + "_" + str(min_len) + "-" + str(max_len) + ".pdf")


# input files and parameters
var_ids_file = sys.argv[1]
intervals = sys.argv[2]
samples_file = sys.argv[3]
outprefix = sys.argv[4]
missed_calls = sys.argv[5]

# interval statistics
interval_to_samples = {}
interval_to_changes = defaultdict(lambda: [])
interval_to_missing = defaultdict(lambda: [])
interval_to_stats = defaultdict(lambda: [0,0])

# all samples to be considered
all_samples = [s.strip() for s in open(samples_file, 'r')]

# alleles that made it into filtered set
filtered_set = {}

# parse the intervals
for line in open(intervals, 'r'):
	fields = line.split('\t')
	samples = [s.strip() for s in fields[9].split(',')]
	interval_to_samples[fields[3]] = samples
	interval_to_changes[fields[3]] = [[0,0] for i in all_samples]
	interval_to_missing[fields[3]] = [[0,0] for i in all_samples]

sys.stderr.write("Read " + str(len(interval_to_samples)) + " gene-interruptive SVs from input file.\n")

# parse filters
for line in open(var_ids_file, 'r'):
	if "variant_id" in line:
		header = line.strip().split()
		continue
	fields = line.split()
	index_id = header.index('variant_id')
	index_conf = header.index('confidence_level')
	id = fields[index_id]
	conf = int(fields[index_conf])
	if conf >= 1:
		filtered_set[id] = True

sys.stderr.write('Read ' + str(len(filtered_set)) + ' filtered IDs from ' + var_ids_file + '.\n')

sample_to_index = {}
for line in sys.stdin:
	if line.startswith('##'):
		continue
	fields = line.strip().split()
	if line.startswith('#'):
		for i,s in enumerate(fields[9:]):
			sample_to_index[s] = i + 9
		continue
	length_change = len(fields[4]) - len(fields[3])	
	info_fields = {v.split('=')[0] : v.split('=')[1] for v in fields[7].split(';') if '=' in v}
	variant_id = info_fields['ID'] if 'ID' in info_fields else fields[2]

	if 'INTERVAL' in info_fields:
		for interval in info_fields['INTERVAL'].split(','):
			interval_to_stats[interval][0] += 1
			if variant_id in filtered_set:
				interval_to_stats[interval][1] += 1

				for i,sample in enumerate(all_samples):
					if not sample in sample_to_index:
						sys.stderr.write("Skip sample " + sample + " since it is not present in VCF.\n")
						continue
					alleles = fields[sample_to_index[sample]].split('|')
					assert len(alleles) == 2
					if alleles[0] == '1':
						sys.stderr.write('\t'.join([interval, sample, '0', variant_id, str(length_change)]) + '\n')
						interval_to_changes[interval][i][0] += length_change
					if alleles[1] == '1':
						sys.stderr.write( '\t'.join([interval, sample, '1', variant_id, str(length_change)]) + '\n' )
						interval_to_changes[interval][i][1] += length_change
					if alleles[0] == '.':
						interval_to_missing[interval][i][0] += 1
					if alleles[1] == '.':
						interval_to_missing[interval][i][1] += 1

sample_header = []
for s in all_samples:
	sample_header.append(s)
header = ['interval_id', 'nr_variants_MC-panel', 'nr_variants_genotyped', 'variants_genotyped[%]'] + sample_header +  ['status', 'mc_filtered_out']
print('\t'.join(header))


missed_mc = defaultdict(lambda: False)
for line in open(missed_calls, 'r'):
	m_id = line.strip()
	missed_mc[m_id] = True

interval_to_status = {}
for interval_id,interval in enumerate(interval_to_samples):
	percentage = (interval_to_stats[interval][1] / max(interval_to_stats[interval][0],1)) * 100.0
	sample_stats = []
	var_present = False
	status = 'NO_MATCHING_CALLS'

	for i in range(len(all_samples)):
		m_h1 =  (interval_to_missing[interval][i][0] / interval_to_stats[interval][1]) if interval_to_stats[interval][1] > 0 else 1.0
		m_h2 =  (interval_to_missing[interval][i][1] / interval_to_stats[interval][1]) if interval_to_stats[interval][1] > 0 else 1.0
		sample = all_samples[i]

		mark = ""
		to_consider = False
		if sample in interval_to_samples[interval]:
			to_consider = True
			mark = "(C)"
		sample_stats.append(str(interval_to_changes[interval][i][0]) + ',' + str(interval_to_changes[interval][i][1]) + ";" + str(m_h1) + ',' + str(m_h2) + mark )
		interval_length = int(interval.split('_')[-1])
		if not to_consider:
			continue
		for change in [interval_to_changes[interval][i][0], interval_to_changes[interval][i][1]]:
			if is_fully_called(interval_length, change):
				# at least one sample has the variant
				status = 'FULLY_CALLED'
			elif status == 'NO_MATCHING_CALLS':
				if is_partially_called(interval_length, change):
					status = 'PARTIALLY_CALLED'

	if not interval_to_stats[interval][0] > 0:
		assert status == "NO_MATCHING_CALLS"
		status = "NO_CALLS_IN_VCF"

	interval_to_status[interval] = status
	to_print = [interval, str(interval_to_stats[interval][0]), str(interval_to_stats[interval][1]), str(percentage), '\t'.join(sample_stats), status, str(missed_mc[interval])]
	print('\t'.join(to_print))

plot(interval_to_samples, interval_to_changes, interval_to_status, 0, 500, outprefix, missed_mc)
plot(interval_to_samples, interval_to_changes, interval_to_status, 500, 5000, outprefix, missed_mc)
plot(interval_to_samples, interval_to_changes, interval_to_status, 5000, 50000, outprefix, missed_mc)
plot(interval_to_samples, interval_to_changes, interval_to_status, 50000, 500000, outprefix, missed_mc)
plot(interval_to_samples, interval_to_changes, interval_to_status, 500000, 5000000, outprefix, missed_mc)
