import sys
from variantclassifier import VariantType, determine_variant_from_line

import argparse

parser = argparse.ArgumentParser(prog='extract-varianttype.py', description=__doc__)
parser.add_argument('vartype', metavar='TYPE', help='variant type. Possible options: snp|small-insertion|small-deletion|midsize-insertion|midsize-deletion|large-insertion|large-deletion|indel|sv|small|midsize|large|all')
args = parser.parse_args()

string_to_type = {
	"snp": [VariantType.snp],
	"small-insertion": [VariantType.small_insertion],
	"small-deletion": [VariantType.small_deletion],
	"small-complex": [VariantType.small_complex],
	"midsize-insertion": [VariantType.midsize_insertion],
	"midsize-deletion": [VariantType.midsize_deletion],
	"midsize-complex": [VariantType.midsize_complex],
	"large-insertion": [VariantType.large_insertion],
	"large-deletion": [VariantType.large_deletion],
	"large-complex": [VariantType.large_complex],
	"indel": [VariantType.snp, VariantType.small_insertion, VariantType.small_deletion, VariantType.small_complex, VariantType.midsize_insertion, VariantType.midsize_deletion, VariantType.midsize_complex],
	"indels": [VariantType.small_insertion, VariantType.small_deletion, VariantType.small_complex, VariantType.midsize_insertion, VariantType.midsize_deletion, VariantType.midsize_complex],
	"sv": [VariantType.large_insertion, VariantType.large_deletion, VariantType.large_complex],
	"small": [VariantType.small_insertion, VariantType.small_deletion, VariantType.small_complex],
	"midsize" : [VariantType.midsize_insertion, VariantType.midsize_deletion, VariantType.midsize_complex],
	"large": [VariantType.large_insertion, VariantType.large_deletion, VariantType.large_complex],
	"all" : [VariantType.snp, VariantType.small_insertion, VariantType.small_deletion, VariantType.small_complex, VariantType.midsize_insertion, VariantType.midsize_deletion, VariantType.midsize_complex, VariantType.large_insertion, VariantType.large_deletion, VariantType.large_complex]
}

assert args.vartype in string_to_type
type_to_extract = string_to_type[args.vartype]

total = 0
written = 0

for line in sys.stdin:
	if line.startswith('#'):
		print(line.strip())
		continue
	vartype = determine_variant_from_line(line)
	fields = line.strip().split()
	is_symbolic_vcf = fields[4].startswith('<')
	if vartype in type_to_extract or is_symbolic_vcf:
		print(line.strip())
		written += 1
	total += 1
sys.stderr.write('Extracted ' + str(written) + ' of ' + str(total) + ' variants.\n')
