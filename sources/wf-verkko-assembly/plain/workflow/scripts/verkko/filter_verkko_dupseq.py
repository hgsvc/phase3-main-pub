#!/usr/bin/env python3

import argparse as argp
import contextlib as ctl
import pathlib as pl
import sys

import pandas as pd
import dnaio


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--fasta",
        "--input",
        "-i",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="input",
        help="Path to input file."
    )

    parser.add_argument(
        "--dedup-report",
        "--report",
        "-r",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="report",
        help="Full path to report file (TSV)."
    )

    parser.add_argument(
        "--output",
        "-o",
        default="stdout",
        dest="output"
    )

    parser.add_argument(
        "--verbose",
        "-vv",
        action="store_true",
        default=False,
        help="Print processing summary to stderr. Default: False"
    )

    args = parser.parse_args()

    return args


def main():

    args = parse_command_line()

    report = pd.read_csv(
        args.report, sep="\t", header=0,
        comment="#"
    )

    input_filename = args.input.name
    select_source = report["seq_source"] == input_filename
    select_discard = report["action"] == "discard"

    discard_sequences = set(
        report.loc[select_source & select_discard, "seq_name"].values
    )

    num_discard = len(discard_sequences)

    if args.output in ["stdout", "-", "/dev/stdout", ""]:
        outfile = sys.stdout.buffer
    else:
        outfile = pl.Path(args.output).resolve()
        outfile.parent.mkdir(exist_ok=True, parents=True)

    discarded = 0
    recnum = 0  # needed for case: empty input
    written_records = 0
    with ctl.ExitStack() as exs:

        read_input = exs.enter_context(
            dnaio.open(args.input, mode="r", fileformat="fasta")
        )

        write_output = exs.enter_context(
            dnaio.open(outfile, mode="w", fileformat="fasta")
        )

        for recnum, record in enumerate(read_input, start=1):
            if record.name in discard_sequences:
                discarded += 1
                continue
            written_records += 1
            write_output.write(record)

    if args.verbose:
        report = (
            f"\n\n=== filter_verkko_dupseq report ==="
            f"\nProcessed file: {input_filename}"
            f"\nProcessed records: {recnum}"
            f"\nDiscarded records: {discarded}"
            f"\nWritten records: {written_records}"
            f"\nFlagged 'discard' in report: {num_discard}"
        )
        sys.stderr.write(report)

    assert num_discard == discarded

    return 0


if __name__ == "__main__":
    main()
