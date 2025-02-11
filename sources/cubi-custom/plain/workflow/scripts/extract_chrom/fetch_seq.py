#!/usr/bin/env python3

import argparse as argp
import functools as fnt
import io
import pathlib as pl

import dnaio
import pandas as pd


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--selected-contigs", "-s",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="selected_contigs",
        help="Table of selected contigs (TSV)"
    )

    parser.add_argument(
        "--fasta-files", "-f",
        type=lambda x: pl.Path(x).resolve(strict=True),
        nargs="+",
        dest="fasta_files",
        help="FASTA sequence files."
    )

    parser.add_argument(
        "--output", "-o",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="output",
        help="Output FASTA file."
    )

    args = parser.parse_args()

    return args


def reverse_complement_sequence(rc_table, sequence):
    return sequence[::-1].translate(rc_table)


def read_selected_contigs_file(table_file):

    assert table_file.suffix == ".tsv"
    table_buffer = io.StringIO()
    with open(table_file, "r") as table:
        for line in table:
            if line.startswith("#"):
                continue
            table_buffer.write(line)
    table_buffer.seek(0)
    df = pd.read_csv(table_buffer, sep="\t", header=0)
    # the following: it is possible that a tig is identified solely
    # based on motif hits and has no alignments to chrY (although that
    # is a rare case)
    all_tigs = set(df["tig_name"].values)
    df = df.loc[df["info_type"] == "aln", :].copy()

    tig_infos = {}
    for tig, alns in df.groupby("tig_name")["num_bp_or_hits"]:
        get_idx = alns.idxmax()
        tig_orient = df.loc[get_idx, "orient_or_specific"]
        assert tig_orient in [-1, 1]
        assert tig not in tig_infos
        tig_infos[tig] = tig_orient

    for tig in all_tigs:
        if tig in tig_infos:
            continue
        # see above: no alignment
        # assign default orientation of forward
        tig_infos[tig] = 1

    return tig_infos


def main():

    args = parse_command_line()

    tig_info = read_selected_contigs_file(args.selected_contigs)

    rc_table = str.maketrans(
        {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    )

    revcomp = fnt.partial(reverse_complement_sequence, rc_table)

    tig_seqs = []
    for fasta_file in args.fasta_files:
        with dnaio.open(fasta_file) as fasta:
            for record in fasta:
                if record.name not in tig_info:
                    continue
                orient = tig_info[record.name]
                if orient < 0:
                    tig_seq = revcomp(record.sequence)
                    tig_name = f"{record.name}.rev"
                else:
                    tig_seq = record.sequence
                    tig_name = f"{record.name}.frw"
                tig_seqs.append((tig_name, tig_seq))

    tig_seqs = sorted(tig_seqs, key=lambda t: len(t[1]), reverse=True)

    args.output.parent.mkdir(exist_ok=True, parents=True)
    with dnaio.FastaWriter(args.output) as writer:
        for tig_name, tig_seq in tig_seqs:
            writer.write(tig_name, tig_seq)

    return 0


if __name__ == "__main__":
    main()
