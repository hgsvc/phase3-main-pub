#!/usr/bin/env python3

import argparse as argp
import pathlib as pl
import sys

import pandas as pd
import xopen


def parse_command_line():

    parser = argp.ArgumentParser()
    parser.add_argument(
        "--input",
        "-i",
        "--paf",
        dest="input",
        type=lambda x: pl.Path(x).resolve(strict=True),
        help="Path to PAF contig-to-reference alignment file",
    )

    parser.add_argument(
        "--output",
        "-o",
        "--tsv",
        dest="output",
        type=lambda x: pl.Path(x).resolve(strict=False),
        help="Path to normalized TSV (will be gzipped) output table."
    )

    args = parser.parse_args()

    return args


def get_fixed_column_definitions():

    # https://lh3.github.io/minimap2/minimap2.html
    PAF_FIXED_COLUMN_NAMES = [
        "query_name", "query_length", "query_start", "query_end",
        "align_orient",
        "target_name", "target_length", "target_start", "target_end",
        "align_matching", "align_total", "mapq"
    ]
    PAF_FIXED_COLUMN_TYPES = [
        str, int, int, int, int, str, int, int, int, int, int, int
    ]
    PAF_FIXED_COLUMN_REPLACEMENTS = [
        None, None, None, None,
        {"+": 1, "-": -1, "*": 0},
        None, None, None, None,
        None, None, None
    ]

    PAF_FIXED_COLUMN_MISSING = [
        None, None, -1, -1,
        0,
        "unaligned", -1, -1, -1,
        0, 0, -1
    ]

    iter_def = (
        PAF_FIXED_COLUMN_NAMES,
        PAF_FIXED_COLUMN_TYPES,
        PAF_FIXED_COLUMN_REPLACEMENTS,
        PAF_FIXED_COLUMN_MISSING
    )

    paf_column_defs = {}
    for cn, (name, dtype, replace, missing) in enumerate(zip(*iter_def), start=1):
        paf_column_defs[cn] = {
            "name": name,
            "dtype": dtype,
            "replace": replace,
            "missing": missing
        }

    return paf_column_defs


def get_variable_column_definitions():
    # https://lh3.github.io/minimap2/minimap2.html

    PAF_VARIABLE_COLUMN_KEYS = [
        "tp", "cm", "s1", "s2", "nm", "md",
        "as", "sa", "ms", "nn", "ts", "cg",
        "cs", "dv", "de", "rl",
        "zd"
    ]
    PAF_VARIABLE_COLUMN_NAMES = [
        "tp_align_type", "cm_num_chain_miniz",
        "s1_chaining_score", "s2_chaining_score",
        "nm_mismatch_gaps", "md_tag", "as_dp_align_score",
        "sa_suppl_align", "ms_dp_align_segment_max_score",
        "nn_ambig_bases", "ts_transcript_strand",
        "cg_cigar", "cs_tag_diffs", "dv_seq_div",
        "de_gapcomp_seq_div", "rl_len_rep_seed_query",
        "zd_undefined"
    ]
    PAF_VARIABLE_COLUMN_TYPES = [
        int, int,
        int, int,
        int, str, int,
        str, int,
        int, int,
        str, str, float,
        float, int,
        int
    ]
    PAF_VARIABLE_COLUMN_REPLACEMENTS = [
        {"P": 1, "p": 1, "S": 2, "s": 2, "I": -1, "i": -1}, None,
        None, None,
        None, None, None,
        None, None,
        None, {"+": 1, "-": -1, "*": 0},
        None, None, None,
        None, None,
        None
    ]

    PAF_VARIABLE_COLUMN_MISSING = [
        0, -1,
        -1, -1,
        -1, "unaligned", -1,
        "unaligned", -1,
        -1, 0,
        "unaligned", "unaligned", -1,
        -1, -1,
        -1
    ]

    iter_def = (
        PAF_VARIABLE_COLUMN_KEYS,
        PAF_VARIABLE_COLUMN_NAMES,
        PAF_VARIABLE_COLUMN_TYPES,
        PAF_VARIABLE_COLUMN_REPLACEMENTS,
        PAF_VARIABLE_COLUMN_MISSING
    )
    paf_column_defs = {}
    for cn, (key, name, dtype, replace, missing) in enumerate(zip(*iter_def), start=13):
        paf_column_defs[key] = {
            "_num": cn,
            "name": name,
            "dtype": dtype,
            "replace": replace,
            "missing": missing
        }

    return paf_column_defs


def get_paf_column_defs():

    paf_fixed_column_defs = get_fixed_column_definitions()
    paf_variable_column_defs = get_variable_column_definitions()

    merged_defs = dict()
    merged_defs.update(paf_fixed_column_defs)
    merged_defs.update(paf_variable_column_defs)
    assert len(merged_defs) == len(paf_fixed_column_defs) + len(paf_variable_column_defs)

    return merged_defs



def read_alignment_file(file_path):

    paf_column_defs = get_paf_column_defs()
    paf_records = []
    processed_lines = 0
    with xopen.xopen(file_path) as paf_aln:
        for line in paf_aln:
            paf_row = dict()
            for cn, content in enumerate(line.strip().split(), start=1):
                if cn < 13:
                    key = cn
                    value = content
                else:
                    try:
                        key, _, value = content.split(":", 2)
                        key = key.lower()
                    except ValueError:
                        print(content)
                        raise
                column_format = paf_column_defs[key]
                if column_format["replace"] is not None:
                    value = column_format["replace"][value]
                paf_row[column_format["name"]] = column_format["dtype"](value)

            paf_records.append(paf_row)
            processed_lines += 1

    paf = pd.DataFrame.from_records(paf_records)
    assert paf.shape[0] == processed_lines

    for column in paf_column_defs.values():
        col_name = column["name"]
        if col_name not in paf.columns:
            continue
        if pd.isnull(paf[col_name]).any():
            paf[col_name].fillna(column["missing"], inplace=True)
            paf[col_name] = paf[col_name].astype(column["dtype"])

    return paf


def main():

    args = parse_command_line()
    alignments = read_alignment_file(args.input)
    alignments.sort_values(
        ["query_name", "align_total", "align_matching"],
        ascending=[True, False, False],
        inplace=True
    )

    args.output.parent.mkdir(exist_ok=True, parents=True)
    alignments.to_csv(args.output, sep="\t", header=True, index=False)

    return 0


if __name__ == "__main__":
    sys.exit(main())
