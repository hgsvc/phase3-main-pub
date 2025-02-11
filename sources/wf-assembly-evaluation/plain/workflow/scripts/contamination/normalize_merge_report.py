#!/usr/bin/env python3

import argparse as argp
import csv
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
        dest="fasta",
        type=lambda x: pl.Path(x).resolve(strict=True),
        help="Path to FCS-processed FASTA file."
    )
    parser.add_argument(
        "--fcs-report",
        "--report",
        "-r",
        dest="report",
        type=lambda x: pl.Path(x).resolve(strict=True)
    )
    mut = parser.add_mutually_exclusive_group()
    mut.add_argument(
        "--adaptor",
        "-a",
        action="store_true",
        default=False,
        dest="adaptor_report"
    )
    mut.add_argument(
        "--contamination",
        "--gx-contam",
        "-c",
        action="store_true",
        default=False,
        dest="contam_report"
    )
    parser.add_argument(
        "--table",
        "--out-table",
        "-t", "-ot",
        dest="table",
        default=pl.Path(".").joinpath("report.norm.tsv"),
        type=lambda x: pl.Path(x).resolve(),
        help="Path to normalized report in TSV format."
    )
    parser.add_argument(
        "--statistics",
        "--out-stats",
        "-s", "-os",
        dest="statistics",
        default=pl.Path(".").joinpath("report.stats.tsv"),
        type=lambda x: pl.Path(x).resolve(),
        help="Path to statistics summary in TSV format."
    )

    args = parser.parse_args()
    if any([args.adaptor_report, args.contam_report]):
        pass
    else:
        raise ValueError(
            "Exactly one of '--adaptor' or '--contamination' must be set."
        )
    return args



def parse_action_range(report_row):

    if pd.isnull(report_row.action_range):
        action_start = 0
        action_end = report_row.seq_length
    else:
        try:
            start, stop = report_row.action_range.split("..")
            action_start = int(start) - 1
            action_end = int(stop)
        except (AttributeError, ValueError, IndexError):
            raise ValueError(f"Cannot parse action range: {report_row}")
    action_name = report_row.action.split("_")[-1]
    assert action_name in ["EXCLUDE", "TRIM"]
    assert action_start < action_end
    assert action_start >= 0
    return action_name, action_start, action_end


def read_adaptor_report(file_path):

    report_header = [
        "name",
        "seq_length",
        "action",
        "action_range",
        "adaptor_name"
    ]
    final_sort_order = [
        "name", "seq_length", "action",
        "action_start", "action_end",
        "adaptor_name"
    ]

    # TODO --- this fails if '#' is used as separator
    # in any of the identifiers in the input
    df = pd.read_csv(
        file_path, sep="\t", comment="#",
        header=None, names=report_header
    )
    if df.empty:
        df = pd.DataFrame(
            columns=final_sort_order
        )
    else:
        # this seems to be an uninformative, fixed
        # component of the identified contaminants
        df["adaptor_name"] = df["adaptor_name"].str.replace("CONTAMINATION_SOURCE_TYPE_", "")
        df["adaptor_name"] = '"' + df["adaptor_name"] + '"'
        norm_actions = df.apply(parse_action_range, axis=1)
        df.drop(["action", "action_range"], axis=1, inplace=True)

        norm_actions = pd.DataFrame.from_records(
            norm_actions, index=df.index,
            columns=["action", "action_start", "action_end"]
        )
        # default: join index-on-index
        df = df.join(norm_actions)
    df = df[final_sort_order]
    return df


def add_sequences_to_adaptor_report(report, fasta_file):

    known_sequences = dict((row.name, row.seq_length) for row in report.itertuples())

    new_records = []
    with dnaio.open(fasta_file) as fasta:
        for seq_record in fasta:
            try:
                seq_length = known_sequences[seq_record.name]
                assert seq_length == len(seq_record.sequence)
            except KeyError:
                seq_length = len(seq_record.sequence)
                record_default = {
                    "name": seq_record.name,
                    "seq_length": seq_length,
                    "action": "PASS",
                    "action_start": 0,
                    "action_end": seq_length,
                    "adaptor_name": "no-adaptor"
                }
                new_records.append(record_default)
    df = pd.DataFrame.from_records(new_records)
    report = pd.concat([report, df], axis=0, ignore_index=False)
    report.sort_values("name", inplace=True)
    return report


def read_contamination_report(file_path):

    report_header = [
        "name",
        "action_start",
        "action_end",
        "seq_length",
        "action",
        "tax_division",
        "pct_cov_in_range",
        "top_tax_name"
    ]
    final_sort_order = [
        "name", "seq_length", "action",
        "action_start", "action_end",
        "pct_cov_in_range", "top_tax_name",
        "tax_division"
    ]

    # TODO --- this fails if '#' is used as separator
    # in any of the identifiers in the input
    df = pd.read_csv(
        file_path, sep="\t", comment="#",
        header=None, names=report_header
    )
    if df.empty:
        df = pd.DataFrame(
            columns=final_sort_order
        )
    else:
        df["tax_division"] = '"' + df["tax_division"] + '"'
        df["top_tax_name"] = '"' + df["top_tax_name"] + '"'
        df["action_start"] = df["action_start"].astype(int)
        df["action_end"] = df["action_end"].astype(int)
        df["action_start"] -= 1
    df = df[final_sort_order]
    return df


def add_sequences_to_contamination_report(report, fasta_file):

    known_sequences = dict((row.name, row.seq_length) for row in report.itertuples())

    new_records = []
    with dnaio.open(fasta_file) as fasta:
        for seq_record in fasta:
            try:
                seq_length = known_sequences[seq_record.name]
                assert seq_length == len(seq_record.sequence)
            except KeyError:
                seq_length = len(seq_record.sequence)
                record_default = {
                    "name": seq_record.name,
                    "seq_length": seq_length,
                    "action": "PASS",
                    "action_start": 0,
                    "action_end": seq_length,
                    "pct_cov_in_range": 100,
                    "top_tax_name": "no-contaminant",
                    "tax_division": "no-contaminant",
                }
                new_records.append(record_default)
    df = pd.DataFrame.from_records(new_records)
    report = pd.concat([report, df], axis=0, ignore_index=False)
    report.sort_values("name", inplace=True)
    return report


def compute_summary_statistics(norm_report):

    stats = []
    total_num_seqs = norm_report["name"].nunique()
    total_seq_length = norm_report["seq_length"].sum()

    stats.append(
        ("TOTAL", "num_sequences", "any", total_num_seqs)
    )
    stats.append(
        ("TOTAL", "seq_length", "any", total_seq_length)
    )

    for action, records in norm_report.groupby("action"):
        action_num_seqs = records["name"].nunique()
        action_total_length = records["seq_length"].sum()
        action_affect_range = int((records["action_end"] - records["action_start"]).sum())
        stats.append(
            (action, "num_sequences", "any", action_num_seqs)
        )
        stats.append(
            (action, "pct_sequences", "any", round(action_num_seqs/total_num_seqs*100, 2))
        )
        stats.append(
            (action, "seq_length", "any", action_total_length)
        )
        stats.append(
            (action, "pct_total_seq_length", "any", round(action_total_length/total_seq_length*100, 2))
        )
        stats.append(
            (action, "affected_length", "any", action_affect_range)
        )
        stats.append(
            (action, "pct_total_affect_length", "any", round(action_affect_range/total_seq_length*100, 2))
        )
        if action != "PASS":
            if "tax_division" in records:
                sub = records.copy()
                for tax_div, sub_records in sub.groupby("tax_division"):
                    tax_div_num_seq = sub_records["name"].nunique()
                    tax_div_total_length = sub_records["seq_length"].sum()
                    tax_div_affect_range = int((sub_records["action_end"] - sub_records["action_start"]).sum())
                    stats.append(
                        (action, "num_sequences", tax_div, tax_div_num_seq)
                    )
                    stats.append(
                        (action, "pct_sequences", tax_div, round(tax_div_num_seq/total_num_seqs*100, 2))
                    )
                    stats.append(
                        (action, "seq_length", tax_div, tax_div_total_length)
                    )
                    stats.append(
                        (action, "pct_total_seq_length", tax_div, round(tax_div_total_length/total_seq_length*100, 2))
                    )
                    stats.append(
                        (action, "affected_length", tax_div, tax_div_affect_range)
                    )
                    stats.append(
                        (action, "pct_total_affect_length", tax_div, round(tax_div_affect_range/total_seq_length*100, 2))
                    )

            elif "adaptor_name" in records:
                sub = records.copy()
                for contam, sub_records in sub.groupby("adaptor_name"):
                    contam_num_seq = sub_records["name"].nunique()
                    contam_total_length = sub_records["seq_length"].sum()
                    contam_affect_range = int((sub_records["action_end"] - sub_records["action_start"]).sum())
                    stats.append(
                        (action, "num_sequences", contam, contam_num_seq)
                    )
                    stats.append(
                        (action, "pct_sequences", contam, round(contam_num_seq/total_num_seqs*100, 2))
                    )
                    stats.append(
                        (action, "seq_length", contam, contam_total_length)
                    )
                    stats.append(
                        (action, "pct_total_seq_length", contam, round(contam_total_length/total_seq_length*100, 2))
                    )
                    stats.append(
                        (action, "affected_length", contam, contam_affect_range)
                    )
                    stats.append(
                        (action, "pct_total_affect_length", contam, round(contam_affect_range/total_seq_length*100, 2))
                    )
            else:
                raise RuntimeError

    stats = pd.DataFrame.from_records(
        stats, columns=["action", "statistic", "group", "value"]
    )
    return stats


def main():

    args = parse_command_line()

    if args.adaptor_report:
        norm_report = read_adaptor_report(args.report)
        norm_report = add_sequences_to_adaptor_report(
            norm_report, args.fasta
        )
        norm_report["screen"] = "ncbi-fcs-adaptor"
    elif args.contam_report:
        norm_report = read_contamination_report(args.report)
        norm_report = add_sequences_to_contamination_report(
            norm_report, args.fasta
        )
        norm_report["screen"] = "ncbi-fcs-gx-contam"
    else:
        raise RuntimeError
    assert norm_report["name"].nunique() == norm_report.shape[0]
    stats = compute_summary_statistics(norm_report)

    args.table.parent.mkdir(exist_ok=True, parents=True)
    norm_report.to_csv(
        args.table, sep="\t", header=True, index=False,
        quoting=csv.QUOTE_NONE
    )

    args.statistics.parent.mkdir(exist_ok=True, parents=True)
    stats.to_csv(
        args.statistics, sep="\t", header=True, index=False,
        quoting=csv.QUOTE_NONE
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
