#!/usr/bin/env python3

import argparse as argp
import pathlib as pl
import re

import dnaio
import pandas as pd


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--fasta-input",
        "-i", "-f",
        type=lambda x: pl.Path(x).resolve(strict=True),
        nargs="+",
        dest="fasta_input",
        help="Path to FASTA input file(s) [space-separated]"
    )

    parser.add_argument(
        "--output",
        "-bed", "-o", "-b",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="output",
        help="Path to output BED-like file (tab-separated)"
    )

    parser.add_argument(
        "--name",
        "-n",
        default="Ngap",
        type=str,
        dest="name",
        help="Name entry for BED-like output. Default: Ngap"
    )

    args = parser.parse_args()
    return args


def main():

    args = parse_command_line()

    ngaps = re.compile("N+")

    bed_ngaps = []
    for fasta_file in args.fasta_input:
        with dnaio.open(fasta_file) as fasta:
            for record in fasta:
                for ngap in ngaps.finditer(record.sequence, re.IGNORECASE):
                    start, end = ngap.span()
                    length = end - start
                    bed_ngaps.append(
                        (record.name, start, end, args.name, length)
                    )

    bed_columns = ["contig", "start", "end", "name", "length"]
    args.output.parent.mkdir(exist_ok=True, parents=True)

    if not bed_ngaps:
        with open(args.output, "w") as bed:
            _ = bed.write("#")
            _ = bed.write("\t".join(bed_columns) + "\n")
    else:
        df = pd.DataFrame.from_records(
            bed_ngaps, columns=bed_columns
        )
        df.sort_values(["contig", "start"], inplace=True)

        with open(args.output, "w") as bed:
            _ = bed.write("#")
            df.to_csv(bed, sep="\t", header=True, index=False)

    return 0


if __name__ == "__main__":
    main()
