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
        "--input",
        "--fasta",
        "-i",
        "-f",
        dest="input_fasta",
        type=lambda x: pl.Path(x).resolve(strict=True),
        help="Path to FCS-processed FASTA file, with contigs flagged with source."
    )

    parser.add_argument(
        "--adapter-table",
        "-at",
        dest="adapter_table",
        type=lambda x: pl.Path(x).resolve(strict=True),
        help="Path to normalized adapter contamination report in TSV format."
    )

    parser.add_argument(
        "--contamination-table",
        "-ct",
        dest="contam_table",
        type=lambda x: pl.Path(x).resolve(strict=True),
        help="Path to normalized foreign contaminants report in TSV format."
    )

    parser.add_argument(
        "--out-pattern",
        "-op",
        dest="out_pattern",
        type=str,
        help="Split input FASTA by tags by adding the keyword SEQTAG in the output file pattern."
    )

    parser.add_argument(
        "--out-contaminants",
        "-oc",
        dest="out_contam",
        type=lambda x: pl.Path(x).resolve(),
        help="Path to FASTA output file for excluded sequences."
    )

    parser.add_argument(
        "--filter-tags",
        "-ft",
        dest="filter_tags",
        type=str,
        default=["hap1", "hap2", "unassigned", "disconnected", "mito", "rdna"],
        nargs="+",
        help="For this set of sequence tags, put contaminated sequences in the '--out-contaminants' file."
    )

    parser.add_argument(
        "--strip-tags",
        "-st",
        action="store_true",
        default=False,
        help="Remove sequence tags in output. Default: False"
    )

    parser.add_argument(
        "--report",
        "-r",
        action="store_true",
        default=False,
        dest="report",
        help="If true, write a brief summary report to stderr. Default: False"
    )

    args = parser.parse_args()
    return args


def read_screening_reports(args):

    # TODO --- this fails if '#' is used as separator
    # in any of the identifiers in the input
    adapter = pd.read_csv(args.adapter_table, sep="\t", comment="#")
    adapter_pass = set(adapter.loc[adapter["action"] == "PASS", "name"].values)

    contam = pd.read_csv(args.contam_table, sep="\t", comment="#")
    contam_pass = set(contam.loc[contam["action"] == "PASS", "name"].values)

    pass_sequences = contam_pass.intersection(adapter_pass)

    return pass_sequences, adapter, contam


def get_splitfile(out_pattern, seqtag):

    splitfile_name = out_pattern.replace("SEQTAG", seqtag)
    splitfile = pl.Path(splitfile_name).resolve()
    splitfile.parent.mkdir(exist_ok=True, parents=True)
    return splitfile


def get_normalized_entity_name(report, name, name_column):

    entity_name = report.loc[report["name"] == name, name_column].values[0]
    entity_name = entity_name.strip('"').replace(" ", ".")
    assert " " not in entity_name
    return entity_name


def trim_contaminated_sequence(report, name, seqtag, lookup_name, sequence, name_column):

    trim_start, trim_end = report.loc[report["name"] == lookup_name, ["action_start", "action_end"]].values[0]
    seq_length = report.loc[report["name"] == lookup_name, "seq_length"].values[0]
    # dropped this check for NA19129 unassigned-0000498.rdna
    # some residual PacBio adapter sequence early in
    # the sequence but not in the beginning
    # assert trim_start == 0 or trim_end == seq_length, f"TRIM start {trim_end} / TRIM end {trim_end} / seq. length {seq_length}"
    add_name = get_normalized_entity_name(report, lookup_name, name_column)

    assert trim_end <= seq_length
    assert trim_start < trim_end

    trim_length = trim_end - trim_start

    if trim_start == 0:
        discard_seq = sequence[:trim_end]
        trimmed_seq = sequence[trim_end:]
    elif trim_end == seq_length:
        discard_seq = sequence[trim_start:]
        trimmed_seq = sequence[:trim_start]
    else:
        discard_seq = sequence[trim_start:trim_end]
        trimmed_seq = sequence[:trim_start] + sequence[trim_end:]

    assert len(discard_seq) == trim_length, f"Discard seq. {len(discard_seq)} vs trim length {trim_length}"

    discard_header = f"{name}|{seqtag}|TRIM|{trim_start}-{trim_end}|{add_name}"

    return trimmed_seq, discard_seq, discard_header


def process_contaminated_sequence(name, seqtag, sequence, adapter_report, contam_report, filter_tags):

    lookup_name = name
    try:
        contam_action = contam_report.loc[contam_report["name"] == lookup_name, "action"].values[0]
    except IndexError:
        lookup_name = f"{name}.{seqtag}"
        contam_action = contam_report.loc[contam_report["name"] == lookup_name, "action"].values[0]
    adapter_action = adapter_report.loc[adapter_report["name"] == lookup_name, "action"].values[0]

    # NB: spelling issues with adaptor vs adapter ...
    name_column_adaptor = "adaptor_name"
    assert name_column_adaptor in adapter_report.columns
    name_column_contam = "tax_division"
    assert name_column_contam in contam_report.columns

    # EXCLUDE flagged sequences are simply removed
    # from the main sequence output files specified
    # by the list of respective tags
    discard = seqtag in filter_tags
    if contam_action != "PASS":
        # contamination report has priority over adapter report
        if contam_action == "EXCLUDE":
            trimmed_seq = None
            trimmed_header = None
            discard_seq = sequence
            discard_reason = get_normalized_entity_name(contam_report, lookup_name, name_column_contam)
            discard_header = f"{name}|{seqtag}|EXCLUDE|contam|{discard_reason}"
        elif contam_action == "TRIM":
            trimmed_seq, discard_seq, discard_header = trim_contaminated_sequence(
                contam_report, name, seqtag, lookup_name, sequence, name_column_contam
            )
            trimmed_header = name
            discard = False
        else:
            raise ValueError(f"Cannot process foreign contamination action: {contam_action} / {name} / {seqtag}")

    elif adapter_action != "PASS":
        # contamination has a PASS, which implies
        # some residual adapter sequence somewhere
        if adapter_action == "EXCLUDE":
            # quite unlikely for adapter contamination
            trimmed_seq = None
            trimmed_header = None
            discard_seq = sequence
            discard_reason = get_normalized_entity_name(adapter_report, lookup_name, name_column_adaptor)
            discard_header = f"{name}|{seqtag}|EXCLUDE|contam|{discard_reason}"
        elif adapter_action == "TRIM":
            trimmed_seq, discard_seq, discard_header = trim_contaminated_sequence(
                adapter_report, name, seqtag, lookup_name, sequence, name_column_adaptor
            )
            trimmed_header = name
            discard = False
        else:
            raise ValueError(f"Cannot process adapter contamination action: {contam_action} / {name} / {seqtag}")
    else:
        raise ValueError(f"Cannot process contaminated sequence record: {name} / {seqtag}")

    return discard, trimmed_header, trimmed_seq, discard_header, discard_seq


def main():

    args = parse_command_line()

    if args.out_pattern in ["stdout", "-", "", "/dev/stdout", "out"]:
        onefile = sys.stdout.buffer
        splitfiles = None
    elif "SEQTAG" not in args.out_pattern:
        onefile = pl.Path(args.out_pattern)
        onefile.parent.mkdir(exist_ok=True, parents=True)
        splitfiles = None
    else:
        onefile = None
        splitfiles = dict()

    args.out_contam.parent.mkdir(exist_ok=True, parents=True)
    contam_out = args.out_contam

    pass_sequences, adapter_report, contam_report = read_screening_reports(args)

    count_records_in = 0
    count_records_out = 0
    count_contam_out = 0
    count_trim_op = 0

    total_read_seq = 0
    total_written_seq = 0
    total_discard_seq = 0
    total_trimmed_seq = 0

    with ctl.ExitStack() as exs:
        if onefile is not None:
            onefile = exs.enter_context(dnaio.open(onefile, mode="w", fileformat="fasta"))

        contam_out = exs.enter_context(dnaio.open(contam_out, mode="w", fileformat="fasta"))

        input_fasta = exs.enter_context(dnaio.open(args.input_fasta, mode="r"))

        for record in input_fasta:
            count_records_in += 1
            total_read_seq += len(record.sequence)
            try:
                stripped_name, seqtag = record.name.rsplit(".", 1)
            except ValueError:
                seqtag = "untagged"
                stripped_name = record.name
            if args.strip_tags:
                out_name = stripped_name
            else:
                out_name = record.name
            if record.name in pass_sequences:
                if onefile is not None:
                    write_out = onefile
                else:
                    try:
                        write_out = splitfiles[seqtag]
                    except KeyError:
                        splitfile = get_splitfile(args.out_pattern, seqtag)
                        write_out = exs.enter_context(dnaio.open(splitfile, mode="w", fileformat="fasta"))
                        splitfiles[seqtag] = write_out
                write_out.write(out_name, record.sequence)
                count_records_out += 1
                total_written_seq += len(record.sequence)
            else:
                discard, trim_header, trim_seq, discard_header, discard_seq = process_contaminated_sequence(
                    out_name, seqtag, record.sequence,
                    adapter_report, contam_report,
                    args.filter_tags
                )
                if discard:
                    # write discard_seq entry to contaminants file
                    contam_out.write(discard_header, discard_seq)
                    count_contam_out += 1
                    total_discard_seq += len(discard_seq)
                else:
                    if onefile is not None:
                        write_out = onefile
                    else:
                        try:
                            write_out = splitfiles[seqtag]
                        except KeyError:
                            splitfile = get_splitfile(args.out_pattern, seqtag)
                            write_out = exs.enter_context(dnaio.open(splitfile, mode="w", fileformat="fasta"))
                            splitfiles[seqtag] = write_out

                    if trim_header is None:
                        # sequence was not trimmed but excluded;
                        # however, it is dumped to the regular
                        # output file because its sequence tag
                        # is not among those to be filtered.
                        write_out.write(out_name, record.sequence)
                        count_records_out += 1
                        total_written_seq += len(record.sequence)
                    else:
                        assert discard_header is not None
                        # sequence was trimmed, so write trimmed
                        # sequence to output file, and discard seq.
                        # to contaminants
                        write_out.write(trim_header, trim_seq)
                        total_written_seq += len(trim_seq)
                        count_records_out += 1

                        # this is just the part of the sequence
                        # that was flagged / trimmed off
                        contam_out.write(discard_header, discard_seq)
                        total_trimmed_seq += len(discard_seq)
                        count_contam_out += 1
                        count_trim_op += 1

    assert count_records_out + count_contam_out >= count_records_in
    assert total_read_seq == total_written_seq + total_trimmed_seq + total_discard_seq

    if args.report:
        summary = (
            f"\n\n=== filter_contaminants report ==="
            f"\nRecords read from input: {count_records_in} / {total_read_seq} bp"
            f"\nRecords written to output: {count_records_out} / {total_written_seq} bp"
            f"\nRecords separated as contaminated: {count_contam_out} / {total_discard_seq} bp"
            f"\nRecords with trimmed sequences: {count_trim_op} / {total_trimmed_seq} bp\n\n"
        )
        sys.stderr.write(summary)

    return 0


if __name__ == "__main__":
    sys.exit(main())
