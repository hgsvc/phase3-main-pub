import sys
from variantclassifier import VariantType, determine_variant_type

region_name = sys.argv[1]

snps = 0
indels = 0
sv_del = 0
sv_ins = 0
sv_complex = 0
all = 0

for line in sys.stdin:
	if line.startswith('#'):
		continue

	fields = line.split()

	# make sure vcf is biallelic
	assert(len(fields[4].split(',')) == 1)	

	info_field = {v.split('=')[0] : v.split('=')[1] for v in fields[7].split(';') if '=' in v}
	assert 'ID' in info_field

	var_id = info_field['ID']
	var_len = int(info_field['ID'].split('-')[-1])
	if 'SNV' in var_id:
		snps += 1
	elif (var_len < 50):
		indels += 1
	elif 'INS' in var_id:
		sv_ins += 1
	elif 'DEL' in var_id:
		sv_del += 1
	else:
		assert 'COMPLEX' in var_id
		sv_complex += 1
	all += 1

print('\t'.join(['region', 'snps', 'indels', 'sv_del', 'sv_ins', 'sv_complex', 'all']))
print('\t'.join([region_name, str(snps), str(indels), str(sv_del), str(sv_ins), str(sv_complex), str(all)]))
