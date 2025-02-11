import sys

def get_type_from_alleles(ref_allele, alt_allele):
	"""
	Determine variant type from REF and ALT allele.
	"""
	if len(ref_allele) < len(alt_allele):
		return "INS"
	elif len(ref_allele) > len(alt_allele):
		return "DEL"
	else:
		return "OTHER"

def get_type_from_ID(varid):
	fields = varid.split('-')
	vartype = fields[2]
	assert vartype in ["SNV", "DEL", "INS", "COMPLEX"]
	return vartype


def get_len_from_ID(varid):
	fields = varid.split('-')
	varlen = fields[-1]
	return varlen	

if __name__ == "__main__":

	# keep track whether header line contains SVTYPE, SVLEN, END tag descriptions
	header_svtype = False
	header_svlen = False
	header_end = False

	for line in sys.stdin:
		if line.startswith('##'):
			print(line.strip())
			if "<ID=SVLEN," in line:
				header_svlen = True
			if "<ID=SVTYPE," in line:
				header_svtype = True
			if "<ID=END," in line:
				header_end = True
			continue
		if line.startswith('#'):
			# add description lines in case they are not present yet
			if not header_svlen:
				print("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of structural variation\">")
			if not header_svtype:
				print("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variation\">")
			if not header_end:
				print("##INFO=<ID=END,Number=1,Type=String,Description=\"End position of structural variation\">")
			print(line.strip())
			continue
		fields = line.strip().split()
		ref_allele = fields[3]
		alt_alleles = fields[4].split(',')
		
		# make sure VCF is biallelic
		if len(alt_alleles) > 1:
			raise Exception("Input VCF is not biallelic. Convert to a biallelic VCF using: bcftools norm -m -any")

		alt_allele = alt_alleles[0]
		is_symbolic = ("<" in alt_allele) and (">" in alt_allele)
		info_fields = {f.split('=')[0] : f.split('=')[1] for f in fields[7].split(';') if '=' in f}
		svtype_present = "SVTYPE" in info_fields
		svlen_present = "SVLEN" in info_fields
		end_present = "END" in info_fields
		
		if is_symbolic:
			# for symbolic variants, make sure that SVTYPE, SVLEN and END tags are present already
			if not end_present:
				raise Exception("Symbolic record at position " + fields[0] + ":" + fields[1] + " misses END tag." )
			if not svtype_present:
				raise Exception("Symbolic record at position " + fields[0] + ":" + fields[1] + " misses SVTYPE tag." )
			if not svlen_present:
				if "DEL" in fields[4]:
					# for DEL, SVLEN can be computed from END and POS fields
					info_fields["SVLEN"] = str(int(info_fields["END"]) - int(fields[1]))
				else:
					raise Exception("Symbolic record at position " + fields[0] + ":" + fields[1] + " misses SVLEN tag." )
		else:
			# explicit vcf record. Add SVTYPE and SVLEN tags in case they are not present already. END tag is not needed,
			# as it can be directly derived from start and ALT allele.
			if not svlen_present:
				sys.stderr.write("Added SVLEN for record " + fields[0] + ":" + fields[1] + "\n")
				assert 'ID' in info_fields
				info_fields["SVLEN"] = get_len_from_ID(info_fields['ID'])
			if not svtype_present:
				sys.stderr.write("Added SVTYPE for record " + fields[0] + ":" + fields[1] + "\n")
				assert 'ID' in info_fields
				info_fields["SVTYPE"] = get_type_from_ID(info_fields['ID'])
		
		fields[7] = ';'.join([k + "=" + v for k,v in info_fields.items()])
		print('\t'.join(fields))
