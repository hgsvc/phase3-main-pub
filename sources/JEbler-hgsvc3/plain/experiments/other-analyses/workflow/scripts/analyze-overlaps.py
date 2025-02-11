import sys
from collections import defaultdict
import gzip

def get_ids_from_vcf(vcf, matches, missed, variant_to_info):
	for line in gzip.open(vcf, 'rt'):
		if line.startswith('#'):
			continue
		fields = line.strip().split()
		info_fields = {v.split('=')[0] : v.split('=')[1] for v in fields[7].split(';') if '=' in v}
		assert 'ID' in info_fields
		var_id = info_fields['ID']
		length = abs(int(info_fields['SVLEN'])) if 'SVLEN' in info_fields else int(var_id.split('-')[-1])
		if length < 50:
			continue
		if not var_id in matches:
			missed.append(var_id)
		end = str(int(fields[1]) + length)
		variant_to_info[var_id] = [fields[0], fields[1], end, str(length)]


def get_ids_from_tsv(file, matches, missed, variant_to_info):
	for line in gzip.open(file, 'rt'):
		if line.startswith('ID'):
			continue
		fields = line.strip().split()
		var_id = fields[0]
		length = abs(int(fields[5]))
		if length < 50:
			continue
		if not var_id in matches:
			missed.append(var_id)
		end = fields[3]
		variant_to_info[var_id] = [fields[1], fields[2], end, str(length)]



callset1_file = sys.argv[1]
callset2_file = sys.argv[2]
callset1_name = sys.argv[3]
callset2_name = sys.argv[4]
unmatched_file = sys.argv[5]

variants = {}

matches = defaultdict(lambda: [])
match_ids1 = set([])
match_ids2 = set([])

missed_ids1 = []
missed_ids2 = []

# first, read all matches from table
for line in sys.stdin:
	if line.startswith('ID'):
		continue
	fields = line.split('\t')
	# only consider SVs
	varlen = int(fields[0].split('-')[-1])
	if varlen < 50:
		continue
	if fields[2] == "A,B":
		# match
		for i in fields[1].split(','):
			matches[fields[0]].append(i)
			match_ids2.add(i)
		match_ids1.add(fields[0])

sys.stderr.write('Read intersections.\n')

# read all IDs from callsets and keep track unmatched ones
if callset1_file.endswith(".vcf.gz"):
	get_ids_from_vcf(callset1_file, match_ids1, missed_ids1, variants)
else:
	get_ids_from_tsv(callset1_file, match_ids1, missed_ids1, variants)
sys.stderr.write('Read ' + callset1_file + '\n')

if callset2_file.endswith(".vcf.gz"):
	get_ids_from_vcf(callset2_file, match_ids2, missed_ids2, variants)
else:
	get_ids_from_tsv(callset2_file, match_ids2, missed_ids2, variants)
sys.stderr.write('Read ' + callset2_file + '\n')


# write the intersection table
header = [ 'ID', 'chromosome', 'start', 'end', 'length', 'var_type', 'quality', 'read_depth', 'allele_depth', 'in_' + callset1_name, 'in_' + callset2_name, 'ID_' + callset1_name, 'ID_' + callset2_name]
print('\t'.join(header))

for m in matches:
	vartype = m.split('-')[2]
	line = [m, variants[m][0], variants[m][1], variants[m][2], variants[m][3], vartype, '.', '.', '.', 'True', 'True', m, ';'.join(matches[m]) ]
	print('\t'.join(line))

with open(unmatched_file, 'w') as outfile:
	for m in missed_ids1:
		vartype = m.split('-')[2]
		line = [m, variants[m][0], variants[m][1], variants[m][2], variants[m][3], vartype, '.', '.', '.', 'True', 'False', m, 'nan']
		print('\t'.join(line))
		outfile.write(m + '\n')
	


for m in missed_ids2:
	vartype = m.split('-')[2]
	line = [m, variants[m][0], variants[m][1], variants[m][2], variants[m][3], vartype, '.', '.', '.', 'False', 'True', 'nan', m]
	print('\t'.join(line))
