import sys

for line in sys.stdin:
	fields = line.strip().split('\t')
	interval_id = '_'.join(fields[:6])
	fields = fields[:3] + [interval_id] + fields[3:]
	print('\t'.join(fields))
	
