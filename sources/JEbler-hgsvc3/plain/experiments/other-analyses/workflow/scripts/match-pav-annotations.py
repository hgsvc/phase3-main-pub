import sys
from collections import defaultdict

pav_annotated = sys.argv[1]
pav_links = sys.argv[2]
mc_pav_matches = sys.argv[3]

pav_name = sys.argv[5]
mc_name = sys.argv[4]

pav_hg38_to_chm13 = defaultdict(lambda: set([]))
for line in open(pav_links, 'r'):
	if line.startswith("A"):
		continue
	fields = line.strip().split()
	pav_hg38_to_chm13[fields[1]].add(fields[0])

pav_chm13_to_mc = defaultdict(lambda: set([]))
index_pav = None
index_mc = None
for line in open(mc_pav_matches, 'r'):
	fields = line.strip().split()
	if line.startswith("ID"):
		index_pav = fields.index(pav_name)
		index_mc = fields.index(mc_name)
		continue
	if (fields[index_pav] != "nan") and (fields[index_mc] != "nan"):
		for pav_id in fields[index_pav].split(';'):
			pav_chm13_to_mc[pav_id].add(fields[index_mc])


for line in open(pav_annotated, 'r'):
	fields = line.strip().split()
	if line.startswith("#"):
		print('\t'.join(fields + ["in_PAV-CHM13", "ID_PAV-CHM13", "in_PanGenie-genotypes", "ID_PanGenie"]))
		continue
	pav_id = fields[3]
	# determine all matched MC calls
	matched_ids = set([])
	pav_matched = set([])
	for i in pav_hg38_to_chm13[pav_id]:
		pav_matched.add(i)
		for j in pav_chm13_to_mc[i]:
			matched_ids.add(j)
	present = "True" if matched_ids else "False"
	present_pav = "True" if pav_matched else "False"
	line = fields + [present_pav, ",".join(pav_matched), present, ",".join(matched_ids)]
	print('\t'.join(line))
