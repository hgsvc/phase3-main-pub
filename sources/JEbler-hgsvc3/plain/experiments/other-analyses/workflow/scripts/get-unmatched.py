import sys
import argparse

parser = argparse.ArgumentParser(prog='intersect-callsets.py', description=__doc__)
parser.add_argument('-c', '--column', required=True, help='column name.')
parser.add_argument('--write-bed', action='store_true', default=False, help="Write BED output.")
args = parser.parse_args()


column_name = args.column
column_id = None
name_id = None
other_ids = []

for line in sys.stdin:
	print_id = True
	fields = line.strip().split("\t")
	if line.startswith("ID"):
		for i,f in enumerate(fields):
			if f == "ID_" + column_name:
				name_id = i
				continue
			if f == "in_" + column_name:
				column_id = i
				continue
			if f.startswith("in_"):
				other_ids.append(i)
		continue
	assert column_id
	if fields[column_id] == "False":
		print_id = False
		continue
	# check if any other column is True
	for i in other_ids:
		if fields[i] == "True":
			print_id = False
			break
	if print_id:
		if args.write_bed:
			print('\t'.join([fields[1], fields[2], fields[3], fields[name_id]]))
		else:
			print(fields[name_id])
