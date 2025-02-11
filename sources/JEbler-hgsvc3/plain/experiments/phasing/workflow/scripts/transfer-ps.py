import sys

for line in sys.stdin:
	if line.startswith("##"):
		print(line.strip())
		continue
	if line.startswith("#"):
		print("##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set\">")
		print(line.strip())
		continue
	fields = line.strip().split()
	info_fields = {f.split("=")[0] : f.split("=")[1] for f in fields[7].split(";") if "=" in f}
	assert "PS" in info_fields
	fields[8] += ":PS"
	for i in range(9, len(fields)):
		fields[i] += ":" + info_fields["PS"]
	print("\t".join(fields))

