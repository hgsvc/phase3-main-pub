import sys
import gzip

vcf1 = sys.argv[1]
vcf2 = sys.argv[2]


# all IDs in the filtered VCF
filtered_ids = {}

# read the filtered VCF and store all IDs present
for line in gzip.open(vcf2, 'rt'):
	if line.startswith('#'):
		continue
	fields = line.strip().split()
	filtered_ids[fields[2]] = True

# read unfiltered vcf and output coordinates of records not in filtered set
for line in gzip.open(vcf1, 'rt'):
	if line.startswith('#'):
		continue
	fields = line.strip().split()
	if fields[2] not in filtered_ids:
		end = str(int(fields[1]) + len(fields[3]))
		print('\t'.join([fields[0], fields[1], end]))
