#!/usr/bin/env python3

import argparse as argp
import collections as col
import io
import functools as ft
import pathlib as pl
import shutil
import sys

import dnaio
import pandas as pd


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--fasta",
        "-f",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="fasta",
        nargs="+",
        help="Path to FASTA file(s).",
        required=True
    )

    parser.add_argument(
        "--cut-table",
        "-t",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="cut_table",
        help="Path to cut table (tsv)",
        required=True
    )

    parser.add_argument(
        "--add-suffix",
        "-r",
        type=str,
        default="roi",
        dest="suffix",
        help="Add this suffix to output file name"
    )

    parser.add_argument(
        "--dump-separately",
        "-s",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="separate",
        default=pl.Path(".").resolve(strict=True),
        help="Folder to dump sequences one per file.",
    )

    parser.add_argument(
        "--clean-out-dir",
        "-c",
        action="store_true",
        dest="clean_out_dir",
        default=False,
        help="Delete out dir at startup"
    )

    parser.add_argument(
        "--dump-merged",
        "-m",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="merged",
        help="Write all sequences into this file.",
        default="stdout"
    )

    parser.add_argument(
        "--dump-stats",
        "-stats",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="stats",
        help="Dump compositional statistics"
    )

    args = parser.parse_args()

    return args


def get_nucleotide_complement_table():

    nuc_comp = {
        "A": "T", "C": "G", "G": "C", "T": "A",
        "a": "t", "c": "g", "g": "c", "t": "a"
    }
    nuc_comp = str.maketrans(nuc_comp)
    return nuc_comp


def reverse_complement_sequence(complement, seq):
    assert isinstance(seq, str)
    return seq[::-1].translate(complement)


def main():

    args = parse_command_line()
    if args.clean_out_dir and args.separate.is_dir():
        shutil.rmtree(args.separate, ignore_errors=False)
    cut_table = pd.read_csv(args.cut_table, sep="\t", header=0)
    cut_table.sort_values(["sample", "source_assembly", "query"], inplace=True)

    collect_all = io.StringIO()
    composition_all = []

    revcomp = ft.partial(reverse_complement_sequence, get_nucleotide_complement_table())

    processed = set()
    for fasta_file in sorted(args.fasta):
        fasta_name = fasta_file.name
        if fasta_name in processed:
            sys.stderr.write(f"\nWARNING: duplicate input - skipping: {fasta_name}\n")
            continue
        processed.add(fasta_name)

        for row in cut_table.itertuples():
            if not fasta_name.startswith(row.source_assembly):
                continue
            cut_begin = row.cut_query_begin
            cut_end = row.cut_query_end
            cut_length = row.cut_length
            query_name = row.query
            query_orientation = row.query_orientation
            sample = row.sample
            src_assembly = row.source_assembly

            with dnaio.open(fasta_file) as fasta:
                for record in fasta:
                    if record.name != query_name:
                        continue
                    if query_orientation < 0:
                        subseq = record.sequence[cut_begin:cut_end]
                        subseq = revcomp(subseq)
                        orient = "REV"
                    else:
                        subseq = record.sequence[cut_begin:cut_end]
                        orient = "FWD"
                    composition = col.Counter(subseq.upper())
                    new_header = (
                        f"{sample}_{query_name}_{args.suffix}"
                        f"_{orient}-{cut_begin}:{cut_end}"
                    )
                    collect_all.write(f">{new_header}\n{subseq}\n")
                    stats = {
                        "sample": sample,
                        "assembly": src_assembly,
                        "contig": query_name,
                        "cut_begin": cut_begin,
                        "cut_end": cut_end,
                        "cut_length": cut_length,
                        "seq_orient": query_orientation
                    }
                    stats.update(composition)
                    composition_all.append(stats)
                    break

            if args.separate is not None:
                outfile = args.separate.joinpath(f"{src_assembly}.{args.suffix}.fasta.gz")
                file_mode = "w"
                outfile.parent.mkdir(exist_ok=True, parents=True)
                if outfile.is_file():
                    file_mode = "a"
                with dnaio.open(outfile, fileformat="fasta", mode=file_mode, compression_level=9) as fasta:
                    fasta.write(new_header, subseq)

    stats = pd.DataFrame.from_records(composition_all)
    if args.stats is not None:
        args.stats.parent.mkdir(exist_ok=True, parents=True)
        stats.to_csv(args.stats, sep="\t", header=True, index=False)

    if args.merged.name in ["stdout", "-", ""]:
        sys.stdout.write(collect_all.getvalue())
    else:
        args.merged.parent.mkdir(exist_ok=True, parents=True)
        with open(args.merged, "w") as dump:
            _ = dump.write(collect_all.getvalue())

    return 0


if __name__ == "__main__":
    main()
