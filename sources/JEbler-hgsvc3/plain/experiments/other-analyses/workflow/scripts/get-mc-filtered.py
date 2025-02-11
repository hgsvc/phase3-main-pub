import sys
import argparse

parser = argparse.ArgumentParser(prog='hwe.py', description=__doc__)
parser.add_argument('min_an', metavar='AN', help="minimum number of paths a variant must be covered.")
parser.add_argument('--chromosomes', metavar='CHROMOSOMES', nargs='+', default=[], help='Only select these chromosomes.')
args = parser.parse_args()

min_an = float(args.min_an)

samples = []

for line in sys.stdin:
	fields = line.split()
	if line.startswith('##'):
		continue
	if line.startswith('#'):
		# keep samples
		samples = fields[9:]
		continue
	chromosome = fields[0].split('.')[-1]
	an_threshold = min_an * 2 * len(samples)
	select_chromosome = (chromosome in args.chromosomes) or (args.chromosomes == [])
	info_fields = {a.split('=')[0] : a.split('=')[1].strip() for a in fields[7].split(';') if '=' in a}
	nr_alt_alleles = len(fields[4].split(','))
	assert 'AN' in info_fields
	assert 'AC' in info_fields
	ac = sum([int(i) for i in info_fields['AC'].split(',')])
	an = int(info_fields['AN'])
	no_n = not 'N' in fields[3] and not 'N' in fields[4]
	no_missing_alt = all(c in 'CAGTcagt,' for c in fields[4])

	printed = (an >= an_threshold) and select_chromosome and no_n and no_missing_alt and (ac > 0)
	if not printed:
		fields[0] = chromosome
		end = str(int(fields[1]) + len(fields[3]))
		print(fields[0], fields[1], end)
