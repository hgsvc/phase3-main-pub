#!/usr/bin/env python3

import argparse as argp
import contextlib as ctl
import operator as op
import pathlib as pl

import pandas as pd
import pysam

def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--norm-paf",
        "-p",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="paf",
        help="Path to normalized PAF alignment table.",
        required=True
    )

    parser.add_argument(
        "--paf-coord-space",
        "-pcs",
        choices=["hpc", "plain"],
        dest="paf_space",
        help="Coordinate space of PAF alignments (hpc or plain [non-hpc]).",
        required=True
    )

    parser.add_argument(
        "--coord-map",
        "-cmap",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="cmap",
        help="Path to coordinate map compatible with PAF alignments (tabix-indexed).",
        required=True
    )

    parser.add_argument(
        "--paf-coord-expand",
        "-c",
        type=str,
        choices=["target", "query"],
        dest="paf_expand",
        help="Expand coordinate space for PAF target or PAF query."
    )

    parser.add_argument(
        "--out-table",
        "-ot",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="out_table",
        help="Path to output table (TSV) with expanded coordinate system."
    )

    args = parser.parse_args()

    return args


def prepare_paf_alignments(paf_file):

    retain_columns = [
        "query_name", "query_length", "query_start", "query_end",
        "target_name", "target_length", "target_start", "target_end",
        "align_orient", "mapq", "tp_align_type", "cg_cigar"
    ]

    # read alignments, assuming that fits in memory
    paf = pd.read_csv(paf_file, sep="\t", header=0, comment="#", usecols=retain_columns)

    # assumption: if the CIGAR string does not contain
    # I or D operations, the interval mapping should be
    # precise, otherwise it is approximate

    cigar_column = [c for c in paf.columns if c.startswith("cg")]
    if len(cigar_column) == 0:
        paf["prec_coord_translation"] = "unqualified"
    elif len(cigar_column) == 1:

        def determine_precision(cigar_string):
            if "D" in cigar_string or "I" in cigar_string:
                return "approximate"
            return "exact"

        cigar_column = cigar_column[0]
        paf["prec_coord_translation"] = paf[cigar_column].apply(determine_precision)
        paf.drop(cigar_column, axis=1, inplace=True)
    else:
        raise ValueError(f"Ambiguous CIGAR string columns: {cigar_column}")

    return paf


def set_expand_target_columns(expand_target):

    expand_name = f"{expand_target}_name"
    expand_start = f"{expand_target}_start"
    expand_end = f"{expand_target}_end"

    return expand_name, expand_start, expand_end


def expand_coordinate_space(paf_row, seq_name, start_coord, end_coord, cmap, get_coords):

    in_start, _, exp_start, _ = get_coords(
        list(cmap.fetch(paf_row[seq_name], paf_row[start_coord], paf_row[start_coord]+1))[0].split()
    )
    _, in_end, _, exp_end = get_coords(
        list(cmap.fetch(paf_row[seq_name], paf_row[end_coord]-1, paf_row[end_coord]))[0].split()
    )
    minus = int(in_start) - paf_row[start_coord]
    plus = int(in_end) - paf_row[end_coord]
    wiggle = f"{minus},{plus}"
    exp_start = int(exp_start)
    exp_end = int(exp_end)
    return exp_start, exp_end, wiggle


def main():

    args = parse_command_line()
    args.out_table.parent.mkdir(exist_ok=True, parents=True)

    paf = prepare_paf_alignments(args.paf)

    exp_name, exp_start, exp_end = set_expand_target_columns(args.paf_expand)
    # if the PAF space is hpc, we expand to plain and vice versa
    exp_space_suffix = "plain" if args.paf_space == "hpc" else "hpc"
    exp_start_column = f"{exp_start}_{exp_space_suffix}"
    exp_end_column = f"{exp_end}_{exp_space_suffix}"

    get_coords = op.itemgetter(
        *tuple([1, 2, 4, 5])
    )
    with ctl.ExitStack() as exs:

        cmap = exs.enter_context(pysam.TabixFile(str(args.cmap)))

        paf[[exp_start_column, exp_end_column, "wiggle"]] = paf.apply(
            expand_coordinate_space,
            axis=1, result_type="expand",
            args=(exp_name, exp_start, exp_end, cmap, get_coords)
        )

        paf.to_csv(args.out_table, header=True, index=False, sep="\t")

    return 0


if __name__ == "__main__":
    main()
