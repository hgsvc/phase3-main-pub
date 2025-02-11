import sys
import matplotlib.pyplot as plt
import numpy as np


len_ins = []
len_del = []
len_snp = []
len_inv = []
sample_name = sys.argv[1]


for line in sys.stdin:
	fields = line.strip().split()
	length = int(fields[3])
	if 'SNV' in fields[4]:
		len_snp.append(1)
	elif 'INS' in fields[4]:
		len_ins.append(length)
	elif 'DEL' in fields[4]:
		len_del.append(length)
	else:
		assert 'INV' in fields[4]
		len_inv.append(length)


plt.figure()

max_v = max(max(len_ins), max(len_del), max(len_inv)) if len_inv else max(max(len_ins), max(len_del))
min_v = 1

for l,k in zip([len_ins, len_del, len_inv, len_snp], ['INS', 'DEL', 'INV', 'SNV']):
	bins = 10 ** np.linspace(np.log10(min_v), np.log10(max_v), 100)
	plt.hist(l, bins = bins, alpha = 0.4, label = k)

plt.title(sample_name)
plt.yscale('symlog')
plt.xscale('symlog')
plt.legend()
plt.xlabel('variant length')
plt.ylabel('Count')
plt.savefig(sys.argv[2])	
