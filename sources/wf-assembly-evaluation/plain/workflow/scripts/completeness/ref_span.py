#!/usr/bin/env python3

import argparse as argp
import io
import pathlib as pl
import sys

import pandas as pd
import xopen


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--approx-align",
        "--input", "-i",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="alignments",
        help="Approximate contig-to-reference alignments (normalized PAF)",
        required=True
    )

    parser.add_argument(
        "--telomere-bed",
        "-telo", "-t",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="telomeres",
        default=None,
        help="Contig telomere annotation produced by seqtk telo [optional]"
    )

    parser.add_argument(
        "--discard-threshold",
        "-dt",
        type=int,
        default=int(5e6),
        dest="discard_threshold",
        help=(
            "Do not discard contigs that span at least this much sequence "
            "on the reference. Default: 5 Mbp"
        )
    )

    parser.add_argument(
        "--ref-span",
        "--output", "-o",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="ref_span",
        help=(
            "Path to output file (BED-like TSV table). "
            "If telomere annotations are provided, includes "
            "indicator columns if contigs are spanning from "
            "telomere/p to telomere/q."
        )
    )

    args = parser.parse_args()

    return args


def load_table_into_buffer(table_file):
    """This function exists for the annoying condition
    that someone used the hash '#' inside of a FASTA
    header as a delimiter, which cannot be handled
    by Pandas' comment policy in read_csv:
    > must be single character, stops parsing anywhere in the line

    Hence, this function loads the entire table into a buffer
    while ignoring only lines that start with a hash #.

    Args:
        table_file (pathlib.Path): File path to alignments table
    """

    table_buffer = io.StringIO()
    with xopen.xopen(table_file, "r") as table:
        for line in table:
            if line.startswith("#"):
                continue
            if not line.strip():
                continue
            table_buffer.write(line)
    # important to reset the buffer to position 0,
    # otherwise pandas.read_csv() will result in
    # empty dataframe
    table_buffer.seek(0)
    return table_buffer


def check_table_is_likely_empty(table):

    num_nan = pd.isnull(table).all(axis=0)
    if num_nan.any():
        raise ValueError("Table is likely empty / malformed")
    return


def select_spanning_alignments(alignments):
    """Select reference/chromosome spanning alignments
    in a greedy fashion:
    1) select largest alignment (column: align_total) as anchor region
    2) extend with left/right margins if applicable
    Args:
        alignments (pandas.DataFrame): Mashmap normalized approx. alignments
    """

    try:
        df = pd.read_csv(alignments, sep="\t", comment="#", header=0)
        check_table_is_likely_empty(df)
    except ValueError:
        err_msg = (
            f"\nError parsing file: {alignments}\n"
            "Standard parsing with pandas.read_csv(... comment='#' ...) failed.\n"
            "Assuming that '#' was used as part of an identifier in the file.\n"
            "Loading table into buffer and restart...\n\n"
        )
        sys.stderr.write(err_msg)
        table_buffer = load_table_into_buffer(alignments)
        df = pd.read_csv(table_buffer, sep="\t", header=0)

    main_chroms = df["target_name"].str.match("[chr0-9A-Z]+$")
    df = df.loc[main_chroms, :].copy()

    span_regions = []
    for qry, aligns in df.groupby("query_name"):

        qry_length = aligns["query_length"].iloc[0]
        trg_length = aligns["target_length"].iloc[0]

        trg_qry_ratio = round(qry_length / trg_length, 2)

        anchor_idx = aligns["align_total"].idxmax()
        anchor_orient = aligns.loc[anchor_idx, "align_orient"]
        anchor_ref_start = aligns.loc[anchor_idx, "target_start"]
        anchor_ref_end = aligns.loc[anchor_idx, "target_end"]
        target_name = aligns.loc[anchor_idx, "target_name"]

        if anchor_orient > 0:
            left_margin = aligns.loc[anchor_idx, "query_start"]
            right_margin = qry_length - aligns.loc[anchor_idx, "query_end"]
        elif anchor_orient < 0:
            left_margin = qry_length - aligns.loc[anchor_idx, "query_end"]
            right_margin = aligns.loc[anchor_idx, "query_start"]
        else:
            raise ValueError(f"No alignment orientation specified: {aligns.loc[anchor_idx, :]}")

        span_from = max(0, anchor_ref_start - left_margin)
        span_to = min(trg_length, anchor_ref_end + right_margin)

        span_regions.append(
            (target_name, span_from, span_to, qry, anchor_orient, trg_qry_ratio, qry_length, trg_length)
        )

    span_regions = pd.DataFrame.from_records(
        span_regions,
        columns=[
            "chrom", "start", "end", "contig",
            "orientation", "length_ratio",
            "contig_length", "chrom_length"
        ]
    )

    return span_regions


def is_contained(check_iv, other_iv):
    """We want to get rid of spurious alignments that are
    fully contained in other alignments to the same region
    (since the list of alignments is sorted by size in
    descending order, subsequent alignments are always
    same size or [usually] smaller).

    Args:
        check_iv (_type_): _description_
        other_iv (_type_): _description_

    Returns:
        bool: _description_
    """

    is_contained = False
    for iv in other_iv:
        if iv.start <= check_iv.start < check_iv.end <= iv.end:
            is_contained = True
            break
    return is_contained


def is_small(iv, lower_size_threshold):
    iv_size = iv.end - iv.start
    return iv_size < lower_size_threshold


def drop_contained_regions(span_regions, discard_threshold):
    """In particular for acrocentric chromosomes,
    the short arms will often be shattered and many
    long-ish sequences cross-align. The greedy selection
    strategy to select the primary "anchor" region for
    each contig is oblivious to more complex alignment
    arrangements and thus, this function simply discards
    all alignments that are fully contained in larger
    ones (overlapping alignments are kept as is).

    Args:
        span_regions (_type_): _description_
    """

    final_regions = []

    for chrom, regions in span_regions.groupby("chrom"):
        if regions.shape[0] == 1:
            final_regions.append(regions.copy())
            continue

        chrom_regions = []

        by_size = regions.sort_values("contig_length", ascending=False, inplace=False)
        for row in by_size.itertuples(index=False):
            if not chrom_regions:
                chrom_regions.append(row)
                continue
            if is_contained(row, chrom_regions) and is_small(row, discard_threshold):
                continue
            chrom_regions.append(row)

        chrom_regions = pd.DataFrame.from_records(
            chrom_regions, columns=span_regions.columns
        )
        final_regions.append(chrom_regions)

    final_regions = pd.concat(final_regions, axis=0, ignore_index=False)

    final_regions.sort_values(["chrom", "start"], inplace=True)

    return final_regions


def read_telomere_annotation(telo_table):

    telo_lut = dict()
    with open(telo_table, "r") as table:
        # format is simple: contig - start - end - size
        for line in table:
            if not line.strip():
                continue
            columns = line.strip().split()
            assert len(columns) == 4
            if int(columns[1]) == 0:  # start
                telo_lut[(columns[0], "begin")] = int(columns[2])
            else:
                assert int(columns[2]) == int(columns[3])
                # end - start = region length
                num_bp = int(columns[2]) - int(columns[1])
                telo_lut[(columns[0], "end")] = num_bp
    return telo_lut


def main():

    args = parse_command_line()

    span_regions = select_spanning_alignments(args.alignments)

    span_regions = drop_contained_regions(span_regions, args.discard_threshold)

    if args.telomeres is not None:
        telo_lut = read_telomere_annotation(args.telomeres)

        # NB: seqtk telo scans simply from left to right, hence
        # the assignment to p- or q-arm can only happen when the
        # alignment orientation is known

        span_regions["telo_p"] = 0
        span_regions["telo_q"] = 0

        in_fwd = span_regions["orientation"] > 0
        span_regions.loc[in_fwd, "telo_p"] = span_regions.loc[in_fwd, "contig"].apply(
            lambda x: telo_lut.get((x, "begin"), 0)
        )
        span_regions.loc[in_fwd, "telo_q"] = span_regions.loc[in_fwd, "contig"].apply(
            lambda x: telo_lut.get((x, "end"), 0)
        )

        in_rev = span_regions["orientation"] < 0
        span_regions.loc[in_rev, "telo_p"] = span_regions.loc[in_rev, "contig"].apply(
            lambda x: telo_lut.get((x, "end"), 0)
        )
        span_regions.loc[in_rev, "telo_q"] = span_regions.loc[in_rev, "contig"].apply(
            lambda x: telo_lut.get((x, "begin"), 0)
        )

    args.ref_span.parent.mkdir(exist_ok=True, parents=True)
    span_regions.to_csv(args.ref_span, sep="\t", header=True, index=False)

    return 0


if __name__ == "__main__":
    main()
