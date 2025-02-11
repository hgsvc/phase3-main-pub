import sys
from collections import defaultdict
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as cm
import argparse

parser = argparse.ArgumentParser(prog='plot-multiway.py', description=__doc__)
parser.add_argument('-tsv', metavar='TSV', nargs='+', required=True, help='Output of whatshap compare multiway comparison.')
parser.add_argument('-trios', metavar='TRIOS', required=True, help="ped file with trio information.")
parser.add_argument('-outname', metavar='OUTNAME', required=True, help="name of the output file.")
parser.add_argument('-title', metavar='TITLE', default="", help="title of the plot.")
args = parser.parse_args()


sample_to_count = defaultdict(lambda: 0)
samples = set([])
sets = set([])

for tsvfile in args.tsv:
	for line in open(tsvfile, 'r'):
		fields = line.strip().split()
		if line.startswith('#'):
			continue
		sample_to_count[(fields[0], fields[2]+fields[3])] += int(fields[4])
		samples.add(fields[0])
		sets.add(fields[2]+fields[3])



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



samples_to_plot = []

unrelated_samples = []
parent_samples = []
child_samples = []

for sample in sorted(samples):
	status = sample_to_status[sample]
	if status == "unrelated":
		unrelated_samples.append(sample)
	elif status == "parent":
		parent_samples.append(sample)
	elif status == "child":
		child_samples.append(sample)

samples_to_plot = unrelated_samples + parent_samples + child_samples
samples_label = [u + " (unrelated)" for u in unrelated_samples] + [p + " (parent)" for p in parent_samples] + [c + " (child)" for c in child_samples]
x_values = [i for i in range(len(samples_to_plot))]

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
fig.subplots_adjust(hspace=0.05)  # adjust space between axes

plt.style.use('tableau-colorblind10')
colors = [cm.to_hex(plt.cm.tab10(i)) for i in range(20)]

set_to_color = {}

sets = sorted(list(sets))

for i,s in enumerate(sets):
	set_to_color[s] = colors[i]


# plot the same data on both axes
for s in sets:
	values = []
	for sample in samples_to_plot:
		values.append(sample_to_count[(sample, s)])
	ax1.plot(x_values, values, label=s, color=set_to_color[s])
	ax2.plot(x_values, values, label=s, color=set_to_color[s])


# zoom-in / limit the view to different portions of the data
ax1.set_ylim(1500000, 3000000)  # outliers only
ax2.set_ylim(0, 60000)  # most of the data

ax1.xaxis.tick_top()
ax1.tick_params(labeltop=False)  # don't put tick labels at the top
ax2.xaxis.tick_bottom()

plt.xticks(x_values, samples_label, rotation=90)
ax1.set_ylabel('Count')
ax2.set_ylabel('Count')
ax1.legend()
ax1.set_title(args.title)
plt.tight_layout()
plt.savefig(args.outname)
