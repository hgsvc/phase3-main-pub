import sys
import pysam

def main():
	bamfile=sys.argv[1]

	bam=pysam.Samfile(bamfile)
	length=0
	chromosomes=['chr' + str(i) for i in range(0,23)] + [str(i) for i in range(0,23)] + ['chrX', 'chrY', 'X', 'Y']
	for e in bam.header.get('SQ'):
		if e['SN'] in chromosomes:
			length += e['LN']
	if length < 1:
		raise ValueError("genome length is 0.")
		sys.exit(1)

	print(length)
	sys.exit(0)
	
main()
