#!/usr/bin/env python3

import argparse as argp
import collections as col
import pathlib as pl

import pandas as pd


####################################
# The following documentation is mostly
# extracted from the paftools.js source
# code and may thus be wrong. Part of it
# is available online here:
# https://lh3.github.io/2020/12/25/evaluating-assembly-quality-with-asmgene
#
# paftools.js output prefixes (lines):
# - full_sgl /// count of full-length single-copy genes
# - full_dup / D / false duplications, should be single-copy
# - frag / F / fragmented genes
# - part50+ / 5 / more than 50% of the transcript aligned
# - part10+ / 1 / more than 10% of the transcript aligned
# - part10- / 0 / less than 10% of the transcript aligned, only traces
# - /// M / missing, no traces of transcript
# - dup_cnt --- count of multi-copy genes
# - dup_sum --- unclear/undefined
#
# Simple QC metrics
# missing single-copy genes / MMC = 1 - (dup_cnt.ASM / dup_cnt.REF)
# missing multi-copy genes / MSC = 1 - (full_sgl.ASM / full_sgl.REF)
#########################################

LINE_PREFIX_MAP = {
    "D": "false_dup",
    "F": "fragmented",
    "5": "incomplete_grt50",
    "1": "incomplete_grt10",
    "0": "incomplete_lst10",
    "M": "missing",
    "d": "multicopy"
}

STATS_ENTRY_MAP = {
    "full_sgl": "single_copy_genes",
    "full_dup": "false_dup_genes",
    "frag": "fragmented",
    "part50+": "incomplete_grt50",
    "part10+": "incomplete_grt10",
    "part10-": "incomplete_lst10",
    "dup_cnt": "multi_copy_genes",
    "dup_sum": "duplicated_sum"  ### unclear what this is
}

BED_CLASSES_MAP = {
    "false_dup": "DUP",
    "fragmented": "FRAG",
    "missing": "MISS",
    "multicopy": "MULTCP"
}

def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--input",
        "-i",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="input_file",
        help="Path to paftools.js asmgene output txt file",
        required=True
    )

    parser.add_argument(
        "--bed-issues",
        "-b",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="bed_issues",
        help="Dump incomplete/missing genes in reference coordinates",
        required=False
    )

    parser.add_argument(
        "--stats-out",
        "-o", "-s",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="statistics",
        help="Dump summary statistics table (.TSV)",
        required=False
    )

    args = parser.parse_args()
    return args


def read_transcript_lines(input_file):

    stats_part = []
    bed_output = []
    label_count = col.Counter()
    with open(input_file, "r") as txt:
        for line in txt:
            if line.startswith("H"):
                stats_part = [line.strip()] + txt.readlines()
                break
            columns = line.strip().split()
            ref_chrom = columns[-3]
            ref_start = int(columns[-2])
            ref_end = int(columns[-1])
            if columns[0] == "d":
                transcript_id = columns[3]
            else:
                transcript_id = columns[2]
            line_label = LINE_PREFIX_MAP[columns[0]]
            label_count[line_label] += 1
            bed_class = BED_CLASSES_MAP.get(line_label, "INCMPL")
            bed_output.append(
                (ref_chrom, ref_start, ref_end, bed_class, transcript_id)
            )
    return stats_part, label_count, bed_output


def parse_asmgene_statistics_lines(stats_part):

    stats_records = []
    for stats_line in stats_part:
        if stats_line.startswith("H"):
            _, metric, ref, query = stats_line.strip().split()
            assert metric == "Metric"
            ref_name = pl.Path(ref).name
            qry_name = pl.Path(query).name
            stats_records.append(
                ("reference_name", ref_name, 1)
            )
            stats_records.append(
                ("query_name", qry_name, 1)
            )
        else:
            assert stats_line.startswith("X")
            _, stat_name, ref_cnt, qry_cnt = stats_line.strip().split()
            stat_label = STATS_ENTRY_MAP[stat_name]
            ref_cnt = int(ref_cnt)
            qry_cnt = int(qry_cnt)
            stats_records.append(
                ("reference", stat_label, ref_cnt)
            )
            stats_records.append(
                ("query", stat_label, qry_cnt)
            )
            if stat_label == "single_copy_genes":
                msc = 100 - round(qry_cnt / ref_cnt * 100, 3)
                stats_records.append(
                    ("both", "miss_sc_pct", msc)
                )
            if stat_label == "multi_copy_genes":
                mmc = 100 - round(qry_cnt / ref_cnt * 100, 3)
                stats_records.append(
                    ("both", "miss_mc_pct", mmc)
                )
    return stats_records


def main():

    args = parse_command_line()

    stats_part, label_count, bed_output = read_transcript_lines(args.input_file)

    stats_records = parse_asmgene_statistics_lines(stats_part)
    [
        stats_records.append(
            ("transcript", key, value)
        ) for key, value in label_count.items()
    ]
    stats_df = pd.DataFrame.from_records(
        stats_records, columns=["entity", "statistic", "value"]
    )

    args.statistics.parent.mkdir(exist_ok=True, parents=True)
    stats_df.to_csv(args.statistics, sep="\t", header=True, index=False)

    bed_df = pd.DataFrame.from_records(
        bed_output, columns=["chrom", "start", "end", "name", "transcript"]
    )
    bed_df.sort_values(["chrom", "start"], inplace=True)

    args.bed_issues.parent.mkdir(exist_ok=True, parents=True)
    with open(args.bed_issues, "w") as bed:
        _ = bed.write("#")
        bed_df.to_csv(bed, sep="\t", header=True, index=False)


    return 0


if __name__ == "__main__":
    main()
