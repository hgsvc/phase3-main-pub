import sys

def parse_info(fields):
	"""
	Parse the info field of a VCF file
	"""
	info_fields = { k.split('=')[0] : k.split('=')[1] for k in fields[7].split(';') if '=' in k }
	return info_fields

print('variant_id\tbubble_id')
for line in sys.stdin:
	if line.startswith('#'):
		continue
	fields = line.strip().split()
	info_fields = parse_info(fields)
	assert 'ID' in info_fields
	ids = set([])
	for allele in info_fields['ID'].split(','):
		for id in allele.split(':'):
			ids.add(id)
	bubble_id = fields[0] + '_' + fields[1]
	for id in ids:
		print(id + '\t' + bubble_id)	
