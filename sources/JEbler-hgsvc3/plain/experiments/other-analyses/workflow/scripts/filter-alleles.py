import sys
from collections import defaultdict

ids = defaultdict(lambda: False)

for line in open(sys.argv[1], 'r'):
	fields = line.strip().split()
	ids[fields[3]] = True


for line in sys.stdin:
	if line.startswith('#'):
		continue
	fields = line.split()
	info_fields = {v.split('=')[0] : v.split('=')[1] for v in fields[7].split(';') if '=' in v}
	assert 'ID' in info_fields
	end = str(int(fields[1]) + len(fields[3]))
	if ids[info_fields['ID']]:
		print('\t'.join([fields[0], fields[1], end]))
