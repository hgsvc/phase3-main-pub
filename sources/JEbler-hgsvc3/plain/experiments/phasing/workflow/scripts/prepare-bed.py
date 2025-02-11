import sys

line_counter = 1

for line in sys.stdin:
	if line.startswith('#'):
		continue
	fields = line.strip().split()
	print('\t'.join([fields[0], fields[1], fields[2], str(line_counter)]))
	line_counter += 1
