import sys
import gzip
from collections import defaultdict

def parse_vcf(vcf, id_in_vcf, index):
	nr_ids = 0

	if vcf.endswith('.gz'):
		for line in gzip.open(vcf, 'rt'):
			if line.startswith('#'):
				continue
			fields = line.split()
			info_field = {k.split('=')[0] : k.split('=')[1] for k in fields[7].split(';') if '=' in k}
			assert 'ID' in info_field
			ids = set([])
			for id_s in info_field['ID'].split(','):
				for id in id_s.split(':'):
					ids.add(id)
			for id in ids:
				id_in_vcf[id][index] = True
				nr_ids += 1
		sys.stderr.write('Read ' + str(nr_ids) + ' IDs from ' + vcf + '.\n')
	else:
		for line in open(vcf, 'r'):
			if line.startswith('#'):
				continue
			fields = line.split()
			info_field = {k.split('=')[0] : k.split('=')[1] for k in fields[7].split(';') if '=' in k}
			assert 'ID' in info_field
			ids = set([])
			for id_s in info_field['ID'].split(','):
				for id in id_s.split(':'):
					ids.add(id)
			for id in ids:
				id_in_vcf[id][index] = True
				nr_ids += 1
		sys.stderr.write('Read ' + str(nr_ids) + ' IDs from ' + vcf + '.\n')




biallelic = sys.argv[1]
multiallelic = sys.argv[2]

id_in_vcf = defaultdict(lambda: [False, False])

parse_vcf(biallelic, id_in_vcf, 0)
parse_vcf(multiallelic, id_in_vcf, 1)

for id, vcfs in id_in_vcf.items():
	if (not vcfs[0]) and vcfs[1]:
		raise RuntimeError('ID in graph but not in callset.')
	if vcfs[0] and not vcfs[1]:
		print(id)
