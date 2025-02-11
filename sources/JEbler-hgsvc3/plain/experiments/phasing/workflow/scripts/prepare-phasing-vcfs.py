import sys
import argparse
import gzip


def check_snp_present(chrom, start, ref, alt, snps):
	if (len(ref) > 1) or len(alt) > 1:
		return False
	position = (chrom, start, alt)
	if position in snps:
		assert snps[position][0] == ref
		return True
	return False

def add_id(old_id, new_id):
	if old_id != ".":
		return old_id + "," + new_id
	else:
		return new_id


if __name__ == "__main__":
	parser = argparse.ArgumentParser(prog='prepare-phasing-vcfs.py', description=__doc__)
	parser.add_argument('phased', metavar='PHASED', help='multisample phased VCF.')
	parser.add_argument('genotypes', metavar='GENOTYPES', help='multisample genotyped VCF.')
	parser.add_argument('rare_snps', metavar='SNPs', help='SNPs to be considered.')
	parser.add_argument('genotypes_output', metavar='GENO_OUT', help='Output VCF with genotypes.')
	parser.add_argument('phasing_output', metavar='PHASE_OUT', help='Output VCF with phasing.')
	parser.add_argument('--geno-source', metavar='GENO_SOURCE', default="genotyping", help="name for the genotyped set.")
	parser.add_argument('--phase-source', metavar='PHASE_SOURCE', default="phasing", help="name for the phased set.")
	args = parser.parse_args()


	# read the rare SNP ids into memory
	# (chrom, start, ALT) -> [REF, None/ID]
	rare_snps = {}
	nr_rare_snps = 0

	for line in open(args.rare_snps, 'r'):
		fields = line.strip().split()
		assert len(fields[2]) == 1
		assert len(fields[3]) == 1
		rare_snps[(fields[0], fields[1], fields[3])] = [fields[2], None]
		nr_rare_snps += 1

	sys.stderr.write("Read " + str(nr_rare_snps) + " from " + args.rare_snps + ".\n")

	# filter the genotypes: ignore variants that are present in the phased set
	nr_genotypes = 0
	nr_present = 0
	with open(args.genotypes_output, 'w') as outfile_genotypes:
		for line in gzip.open(args.genotypes, 'rt'):
			if line.startswith('##'):
				# header line
				# skip INFO lines, since INFO field will be set to '.'
				if ("INFO=" in line) or ("bcftools" in line):
					continue
				outfile_genotypes.write(line.strip() + '\n')
				continue
			if line.startswith('#'):
				# print INFO metaline, describing source of the genotypes
				outfile_genotypes.write("##INFO=<ID=SOURCE,Number=1,Type=String,Description=\"Source of the genotyping information.\">\n")
				outfile_genotypes.write(line.strip() + '\n')
				continue
			fields = line.strip().split()
			if len([a for a in fields[4].split(',')]) > 1:
				raise Exception("Genotyped VCF must be biallelic.")
			snp_present = check_snp_present(fields[0], fields[1], fields[3], fields[4], rare_snps)

			if snp_present:
				# store ID
				if fields[2] != '.':
					position = (fields[0], fields[1], fields[4])
					rare_snps[position][1] = fields[2]
					nr_present += 1
			else:
				# set INFO field to '.'
				fields[7] = 'SOURCE=' + args.geno_source
				outfile_genotypes.write('\t'.join(fields) + '\n')
				nr_genotypes += 1

	sys.stderr.write("Wrote " + str(nr_genotypes) + " genotypes.\n")
	sys.stderr.write("Of those, " + str(nr_present) + " SNPs were present in both VCFs.\n")

	nr_phasing = 0
	with open(args.phasing_output, 'w') as outfile_phasing:
		for line in gzip.open(args.phased, 'rt'):
			if line.startswith('##'):
				# header line
				# skip INFO lines, since INFO field will be set to '.'
				if ("INFO=" in line) or ("bcftools" in line):
					continue
				outfile_phasing.write(line.strip() + '\n')
				continue
			if line.startswith('#'):
				# print INFO metaline, describing source of the genotypes
				outfile_phasing.write("##INFO=<ID=SOURCE,Number=1,Type=String,Description=\"Source of the genotyping information.\">\n")
				outfile_phasing.write(line.strip() + '\n')
				continue
			fields = line.strip().split()
			# only consider SNPs
			if (len(fields[3]) > 1) or (len(fields[4]) > 1):
				continue
			if len([a for a in fields[4].split(',')]) > 1:
				raise Exception("Phasing VCF must be biallelic.")
			snp_present = check_snp_present(fields[0], fields[1], fields[3], fields[4], rare_snps)
			if snp_present:
				position = (fields[0], fields[1], fields[4])
				# if there is an additional ID, add it
				if rare_snps[position][1]:
					fields[2] = add_id(fields[2], rare_snps[position][1])
				# set INFO field to '.'
				fields[7] = 'SOURCE=' + args.phase_source
				outfile_phasing.write('\t'.join(fields) + '\n')
				nr_phasing += 1
	sys.stderr.write("Wrote " + str(nr_phasing) + " phased genotypes.\n")
		
