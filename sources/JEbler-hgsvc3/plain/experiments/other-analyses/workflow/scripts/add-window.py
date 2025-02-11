import sys

for line in sys.stdin:
	fields = line.strip().split('\t')
	# add 10 kb left and right
	fields[1] = str(int(fields[1]) - 10000)
	fields[2] = str(int(fields[2]) + 10000)
	print('\t'.join(fields))	
