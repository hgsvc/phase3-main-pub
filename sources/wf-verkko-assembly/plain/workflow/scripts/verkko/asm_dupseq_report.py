#!/usr/bin/env python3

import argparse as argp
import pathlib as pl
import sys

import pandas as pd


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--input",
        "-i",
        nargs="+",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="input",
        help="Path to input files."
    )

    parser.add_argument(
        "--report",
        "-r",
        type=lambda x: pl.Path(x).resolve(),
        default=pl.Path(".").joinpath("report.tsv.gz"),
        dest="report",
        help="Full path to report file (TSV)."
    )

    parser.add_argument(
        "--summary",
        "-s",
        type=lambda x: pl.Path(x).resolve(),
        default=pl.Path(".").joinpath("summary.tsv"),
        dest="summary",
        help="Full path to summary file (TSV)"
    )

    args = parser.parse_args()

    return args


def compute_record_multiplicites(seq_ids):

    columns = ["seq_name", "seq_hash"]
    new_names = ["name_multiplicity", "sequence_multiplicity"]

    for c, n in zip(columns, new_names):
        mult = seq_ids[c].value_counts()
        mult.rename(n, inplace=True)
        seq_ids = seq_ids.merge(mult, left_on=c, right_index=True, how="outer")

    seq_ids.reset_index(inplace=True, drop=True)

    uniq_name = seq_ids["name_multiplicity"] == 1
    uniq_seq = seq_ids["sequence_multiplicity"] == 1
    seq_ids.loc[uniq_name & uniq_seq, "is_unique_record"] = 1
    seq_ids.loc[uniq_name & uniq_seq, "seq_group"] = "unique"
    seq_ids.loc[uniq_name & uniq_seq, "action"] = "keep"


    return seq_ids


def classify_non_unique_record(seq_ids):

    for row in seq_ids.loc[seq_ids["is_unique_record"] < 1, :].itertuples():
        nmult = row.name_multiplicity
        smult = row.sequence_multiplicity
        source_prio = row.source_priority
        # case 1: remove duplicate in prio 2 source
        if nmult > 1 and smult > 1:
            assert nmult == smult
            seq_group = "duplicate"
            if source_prio > 0:
                action = "keep"
                comment = "dupseq.in.high.priority.source"
            elif source_prio < 0:
                action = "discard"
                comment = "dupseq.in.low.priority.source"
            else:
                raise ValueError(row)
        # case 2: tricky, unique sequence
        # but name used several times;
        # assume that is an error
        elif nmult > 1 and smult == 1:
            seq_group = "error"
            action = "raise"
            comment = "duplicated.name.with.unique.sequence"
        # case 3: several names identify the same
        # sequence; weird, but may happen when
        # assembling highly repetitive parts of the
        # genome?
        elif nmult == 1 and smult > 1:
            seq_group = "duplicate"
            action = "keep"
            comment = "unique.name.with.duplicated.sequence"
        else:
            raise ValueError(f"Cannot process row: {row}")

        seq_ids.loc[row.Index, "seq_group"] = seq_group
        seq_ids.loc[row.Index, "action"] = action
        seq_ids.loc[row.Index, "comment"] = comment

    return seq_ids


def main():

    args = parse_command_line()

    seq_ids = []
    for tsv_file in args.input:
        df = pd.read_csv(
            tsv_file, sep="\t", header=0, comment="#",
            usecols=["seq_name", "seq_source", "seq_hash", "seq_length"]
        )
        seq_ids.append(df)
    seq_ids = pd.concat(seq_ids, axis=0, ignore_index=False)
    seq_ids.reset_index(inplace=True, drop=True)

    seq_ids["source_priority"] = 1
    seq_ids.loc[seq_ids["seq_source"].str.contains("disconnected"), "source_priority"] = -1
    seq_ids["is_unique_record"] = 0
    seq_ids["seq_group"] = "undefined"
    seq_ids["action"] = "undefined"
    seq_ids["comment"] = "no.comment"

    seq_ids = compute_record_multiplicites(seq_ids)
    seq_ids = classify_non_unique_record(seq_ids)
    seq_ids.sort_values(["source_priority", "seq_source", "seq_name", "seq_length"], inplace=True)

    if "error" in seq_ids["seq_group"].values:
        for row in seq_ids.loc[seq_ids["seq_group"] == "error", :].itertuples():
            sys.stderr.write(f"\nERROR --- {row}\n")
        raise ValueError("Assembled sequences contain records classified as error")

    summary = seq_ids.groupby(["seq_source", "seq_group"])["action"].value_counts()
    args.summary.parent.mkdir(exist_ok=True, parents=True)
    summary.to_csv(args.summary, index=True, header=True, sep="\t")

    out_column_order = [
        "seq_name", "seq_length", "seq_source",
        "source_priority",
        "seq_group", "action", "comment",
        "is_unique_record",
        "name_multiplicity", "sequence_multiplicity",
        "seq_hash"
    ]
    assert all(c in out_column_order for c in seq_ids.columns)

    args.report.parent.mkdir(exist_ok=True, parents=True)
    seq_ids[out_column_order].to_csv(args.report, index=False, header=True, sep="\t")

    return 0


if __name__ == "__main__":
    main()
