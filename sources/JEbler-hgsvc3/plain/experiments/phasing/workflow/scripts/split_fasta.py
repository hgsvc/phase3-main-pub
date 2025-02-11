import sys
import argparse


def remove_ns(outname):
	with open(outname, 'w') as outfile:
		total_written = 0
		currently_written = 0
		previous_label = ""
		nr_blocks = 0
		prev_n = False
		prev_break = True

		for line in sys.stdin:
			if line == "":
				continue
			if line.startswith(">"):
				previous_label = line.strip().split()[0]
				nr_blocks = 0
				prev_n = True
				continue
			for c in line:
				if c == '\n':
					continue
				if c in ['n', 'N']:
					prev_n = True
					continue
				elif prev_n:
					if not prev_break:
						outfile.write("\n")
					outfile.write(previous_label + "_" + str(nr_blocks) + "\n" + c)
					total_written += 1
					prev_break = False
					currently_written = 1
					nr_blocks += 1
					prev_n = False
				
				else:
					if currently_written % 70 == 0:
						outfile.write("\n")
						prev_break = True
					outfile.write(c)
					prev_break = False
					currently_written += 1
					total_written += 1
		if not prev_break:
			outfile.write("\n")
	print("Wrote " + str(total_written) + " non-N bases to output file.")


if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='split_fasta.py', description=__doc__)
	parser.add_argument('-o', metavar='FASTA', required=True, help='Name of the output FASTA file.')
	args = parser.parse_args()

	remove_ns(args.o)
