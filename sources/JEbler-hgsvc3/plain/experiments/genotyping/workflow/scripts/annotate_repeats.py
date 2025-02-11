import sys
import argparse
import gzip

def parse_info(fields):
	"""
	Parse the info field of a VCF file
	"""
	info_fields = { k.split('=')[0] : k.split('=')[1] for k in fields[7].split(';') if '=' in k }
	return info_fields


def produce_tsv(names):
	print('variant_id\t' + '\t'.join([n + '_overlaps' for n in names]))

	for line in sys.stdin:
		if line.startswith('#'):
			continue
		fields = line.strip().split()
		info_fields = parse_info(fields)
		# make sure VCF is biallelic
		len(fields[4].split(',')) == 1
		if 'ID' in info_fields:
			var_id = info_fields['ID'].strip()
		else:
			var_id = fields[2]
		n = len(names)
		overlaps = [o for o in fields[-n:]]
		print(var_id + '\t' + '\t'.join(overlaps))


def produce_vcf(names, vcf):
	for line in gzip.open(vcf, 'rt'):
		if line.startswith('##'):
			print(line.strip())
			continue
		if line.startswith('#'):
			print(line.strip() + '\t' + '\t'.join([n + '_overlaps' for n in names]))
			break

	for line in sys.stdin:
		if line.startswith('#'):
			continue
		fields = line.strip().split()
		# make sure VCF is biallelic
		len(fields[4].split(',')) == 1
		n = len(names)
		overlaps = [o for o in fields[-n:]]
		print(line.strip())



if __name__== "__main__":
	parser = argparse.ArgumentParser(prog='annotate_repeats.py', description=__doc__)
	parser.add_argument('-vcf', metavar='VCF', help='VCF that was annotated, needed for header lines.')
	parser.add_argument('-names', metavar='NAMES', required=True, nargs='+', help='names of regions that were used for annotation.')
	parser.add_argument('-format', metavar='FORMAT', required=True, choices=['vcf', 'tsv'], help='output VCF format (with added columns) or TSV table.')
	args = parser.parse_args()

	assert args.format in ['vcf', 'tsv']

	if args.format == 'tsv':
		produce_tsv(args.names)
	else:
		produce_vcf(args.names, args.vcf)

