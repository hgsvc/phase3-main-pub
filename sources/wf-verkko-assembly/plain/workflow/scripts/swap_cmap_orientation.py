#!/usr/bin/env python3

import operator as op
import sys


def main():

    # expected line input format
    # ["#name", "hpc_start", "hpc_end", "seqnum", "plain_start", "plain_end", "hpc_ratio"]
    format_line = sys.stdin.readline().strip().split()

    num_columns = len(format_line)
    column_indices = list(range(num_columns))
    column_indices[1] = 4
    column_indices[2] = 5
    column_indices[4] = 1
    column_indices[5] = 2

    swap = op.itemgetter(*tuple(column_indices))
    new_format_line = "\t".join(swap(format_line))

    sys.stdout.write(new_format_line + "\n")

    for line in sys.stdin:
        out_line = "\t".join(swap(line.strip().split()))
        sys.stdout.write(out_line + "\n")

    return


if __name__ == "__main__":
    main()
