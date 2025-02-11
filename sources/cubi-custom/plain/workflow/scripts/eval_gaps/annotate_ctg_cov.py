#!/usr/bin/env python3

import argparse as argp
import pathlib as pl

import pandas as pd


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--ctg-cov", "-c",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="ctg_cov",
        help="TSV table of assembly contig coverage (10 kbp bins) in chm13 coordinates."
    )

    parser.add_argument(
        "--gaps", "-g",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="gaps",
        help="TSV table of common HPRC haps in chm13 coordinates."
    )

    parser.add_argument(
        "--out-agg", "-o",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="output",
        help="Aggregated output per gap (TSV)."
    )

    parser.add_argument(
        "--select-mapq", "-mq",
        type=str,
        choices=["MQ60", "MQ00"],
        default="MQ60",
        dest="select_mapq",
        help="Select MAPQ threshold. Default: MQ60"
    )

    parser.add_argument(
        "--sample", "-s",
        type=str,
        default=None,
        dest="sample",
        help="Sample name to add to output table."
    )

    parser.add_argument(
        "--assembly-sex", "-x",
        type=str,
        choices=["female", "female,male", "male,female", "unknown"],
        default="unknown",
        dest="assembly_sex",
        help="Assembly sex. Default: unknown"
    )

    args = parser.parse_args()

    return args


def load_binned_contig_coverage(file_path, select_mapq):

    ctg_cov = pd.read_csv(file_path, sep="\t", header=[0,1,2], index_col=[0,1,2])
    keep_columns = []
    asm_units = ["asm-hap1", "asm-hap2", "asm-unassigned"]
    known_asm_units = set()
    for c in ctg_cov.columns:
        keep_unit = c[0] in asm_units
        keep_mq = c[1] == select_mapq
        if keep_unit and keep_mq:
            known_asm_units.add(c[0])
            keep_columns.append(c)
    ctg_cov = ctg_cov[keep_columns].copy()
    ctg_cov.columns = ctg_cov.columns.droplevel(["statistic", "mapq"])
    ctg_cov = ctg_cov.reset_index(drop=False, inplace=False)

    return ctg_cov, sorted(known_asm_units)


def main():

    args = parse_command_line()

    gaps = pd.read_csv(args.gaps, sep="\t", header=0)

    ctg_cov, asm_units = load_binned_contig_coverage(args.ctg_cov, args.select_mapq)

    sample_name = "sample" if args.sample is None else args.sample
    assembly_sex = args.assembly_sex
    if "," in assembly_sex:
        hap1_sex, hap2_sex = assembly_sex.split(",")
    else:
        hap1_sex, hap2_sex = assembly_sex, assembly_sex
    au_sex = []
    for au in asm_units:
        if "hap1" in au or "h1" in au:
            au_sex.append(hap1_sex)
        elif "hap2" in au or "h2" in au:
            au_sex.append(hap2_sex)
        else:
            if assembly_sex == "female":
                au_sex.append("female")
            else:
                au_sex.append("unknown")

    aggregate = []
    for row in gaps.itertuples():
        if assembly_sex == "female" and row.chrom == "chrY":
            continue
        select_chrom = ctg_cov["chrom"] == row.chrom
        select_start = ctg_cov["start"] >= row.win_start  # NB: binned / blunt window ends
        select_end = ctg_cov["end"] <= row.win_end
        selector = select_chrom & select_start & select_end

        subset = ctg_cov.loc[selector, :]
        # this cannot be empty, only the coverage can be zero
        num_windows = subset.shape[0]
        for au, sex in zip(asm_units, au_sex):
            if sex == "female" and row.chrom == "chrY":
                continue
            cov_tnz = (subset[au] > 0).sum()
            cov_one = (subset[au] == 1).sum()
            cov_two = (subset[au] == 2).sum()
            record = (
                row.chrom, row.start, row.end,
                row.name, sample_name, sex,
                au, num_windows, cov_tnz, cov_one, cov_two
            )
            aggregate.append(record)

    agg_columns = [
        "chrom", "start", "end",
        "gap_id", "sample", "asm_sex",
        "asm_unit", "num_windows",
        "cov_tnz", "cov_one", "cov_two"
    ]
    aggregate = pd.DataFrame.from_records(
        aggregate, columns=agg_columns
    )
    aggregate.sort_values(["chrom", "start", "end"])

    args.output.parent.mkdir(exist_ok=True, parents=True)
    aggregate.to_csv(args.output, sep="\t", header=True, index=False)

    return 0


if __name__ == "__main__":
    main()
