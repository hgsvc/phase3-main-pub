import sys
from collections import defaultdict
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt



sample_to_concordance = {}
sample_to_stats = defaultdict(lambda: [None, None, None])


for filename in sys.stdin:
	sample = filename.strip().split("/")[-1].split("_")[-1].split(".")[0]
	for line in open(filename.strip(), 'r'):
		fields = line.strip().split()
		if line.startswith("weighted_concordance"):
			concordance = float(fields[-1])
			sample_to_concordance[sample] = concordance * 100.0
		if len(fields) < 3:
			continue
		if (fields[0] == "False") and (fields[1] == "True"):
			sample_to_stats[sample][0] = int(fields[2])
		elif (fields[0] == "True") and (fields[1] == "False"):
			sample_to_stats[sample][1] = int(fields[2])
		elif (fields[0] == "True") and (fields[1] == "True"):
			sample_to_stats[sample][2] = int(fields[2])

# print table with concordance and variant numbers
header = ["sample", "weighted_genotype_concordance", "variants_in_intersection", "variants_only_in_truthset", "variants_only_in_genotyped_set"]
print("\t".join(header))

for sample in sample_to_concordance.keys():
	line = [sample, str(sample_to_concordance[sample]), str(sample_to_stats[sample][2]), str(sample_to_stats[sample][1]), str(sample_to_stats[sample][0])]
	print("\t".join(line))


## plot genotype concordances

concordances = [sample_to_concordance[s] for s in sample_to_concordance.keys()]

plt.figure()
plt.hist(concordances, bins=20, range=[99,100])
plt.xlabel("weighted genotype concordances")
plt.ylabel("count")
plt.title(sys.argv[2])
plt.savefig(sys.argv[1])
