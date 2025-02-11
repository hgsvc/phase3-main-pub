import sys
import argparse

parser = argparse.ArgumentParser(prog='remove-missing.py', description="Remove all positions for which the given sample misses genotype information.")
parser.add_argument('sample', metavar='SAMPLE', help='which sample to consider.')
args = parser.parse_args()

header = None
sample_of_interest = args.sample
total_variants = 0
skipped_variants = 0


for line in sys.stdin:
	if line.startswith('##'):
		print(line.strip())
		continue
	fields = line.split()
	if line.startswith('#'):
		header = fields
		assert sample_of_interest in header
		print(line.strip())
		continue
	total_variants += 1
	sample_index = header.index(sample_of_interest)
	genotype = fields[sample_index]
	if '.' in genotype:
		skipped_variants += 1
		continue
	print(line.strip())
sys.stderr.write('Skipped ' + str(skipped_variants) + '/' + str(total_variants) + ' due to missing genotype information in sample ' + sample_of_interest + '\n')
