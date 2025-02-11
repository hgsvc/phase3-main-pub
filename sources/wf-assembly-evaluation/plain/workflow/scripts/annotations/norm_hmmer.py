#!/usr/bin/env python3

import argparse as argp
import collections as col
import io
import pathlib as pl
import sys

import dnaio
import pandas as pd
import xopen

# this pertains at least to version 3.4
HMMER_TABLE_COLUMNS = [
    ("target", True, str),
    ("target_accession", False, None),
    ("query", True, str),
    ("query_accession", False, None),
    ("query_hit_start", True, int),
    ("query_hit_end", True, int),
    ("target_hit_start", True, int),
    ("target_hit_end", True, int),
    ("target_env_start", True, int),
    ("target_env_end", True, int),
    ("target_length", True, int),
    ("target_strand", True, str),
    ("evalue", True, str),  # force string here to avoid ugly floating point output
    ("bit_score", True, float),
    ("bias", True, float),
    ("description", False, str)
]

HMMER_TABLE_NAMES = [t[0] for t in HMMER_TABLE_COLUMNS]
HMMER_TABLE_USE_COLS = [t[0] for t in HMMER_TABLE_COLUMNS if t[1]]
HMMER_TABLE_DTYPES = dict((t[0], t[2]) for t in HMMER_TABLE_COLUMNS)

HMMER_BED_COLUMNS = [
    "target",
    "target_hit_start",
    "target_hit_end",
    "query",
    "bit_score",
    "target_strand",
    "hit_hiq",
    "evalue"
]

# this info is coded here in case empty
# output needs to be created
HMMER_AGG_COLUMNS = [
    "target", "query", "target_length",
    "num_hits_total", "num_hits_hiq",
    "num_hits_loq", "pct_hits_hiq",
    "num_bp_hiq", "pct_bp_hiq",
    "mean_pct_len_hiq", "median_pct_len_hiq",
    "mean_score_hiq", "median_score_hiq"
]


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--hmmer-table", "-ht", "-tab",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="hmmer_tblout",
        help="Path to nhmmer tblout output file.",
        required=True
    )

    parser.add_argument(
        "--motif-file", "-mf", "-f",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="motif_file",
        help="Path to motif file (FASTA)."
    )

    parser.add_argument(
        "--score-threshold", "-st",
        type=int,
        dest="score_threshold",
        default=-1,
        help=(
            "Specify bit score threshold to label motif matches as high quality. "
            "Set to -1 to disable (no thresholding). Thresholding operation is 'larger than / > T'. "
            "Default: -1"
        )
    )

    parser.add_argument(
        "--evalue-threshold", "-et",
        type=float,
        dest="evalue_threshold",
        default=-1,
        help=(
            "Specify E-value threshold to label motif matches as high quality. "
            "Set to -1 to disable (no thresholding). Thresholding operation is 'smaller than / < T'. "
            "Default: -1"
        )
    )

    parser.add_argument(
        "--norm-table", "-out", "-o",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="norm_table",
        help="Normalized output table (TSV).",
        required=True,
    )

    parser.add_argument(
        "--bed-all", "-bed",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="bed_all",
        help="Normalized motif matches reduced to BED-like output format.",
        required=True,
    )

    parser.add_argument(
        "--bed-hiq", "-hiq",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="bed_hiq",
        help="Normalized motif matches reduced to high-quality hits dumped to BED-like output format.",
        required=True,
    )

    parser.add_argument(
        "--aggregated", "-agg", "-a",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="aggregated",
        help="Aggregated statistics for all motif matches (TSV).",
        required=True,
    )

    parser.add_argument(
        "--force-empty-output", "-empty",
        action="store_true",
        default=False,
        dest="force_empty",
        help=(
            "If the input is empty (no motif hits), created empty output files "
            "(just header lines). Default: False"
        )
    )

    parser.add_argument(
        "--add-empty-flag", "-flag",
        action="store_true",
        default=False,
        dest="empty_flag",
        help=(
            "If an empty output files is created (see option --force-empty-output), "
            "create a flag file (<FILENAME.EXT>.EMPTY) right next to it to indicate "
            "that the empty output is intentional. Default: False"
        )
    )

    args = parser.parse_args()

    return args


def get_motif_info(fasta_file):

    record_counter = 0
    with dnaio.open(fasta_file) as fasta:
        for record in fasta:
            assert record_counter == 0, f"Multiple motifs in input FASTA are not supported: {fasta_file}"
            query_name = record.name
            query_length = len(record.sequence)
            record_counter += 1
    return query_name, query_length


def load_table_into_buffer(table_file):
    """This function exists for the annoying condition
    that someone used the hash '#' inside of a FASTA
    header as a delimiter, which cannot be handled
    by Pandas' comment policy in read_csv:
    > must be single character, stops parsing anywhere in the line

    Hence, this function loads the entire table into a buffer
    while ignoring only lines that start with a hash #.

    Args:
        table_file (pathlib.Path): File path to HMMER output table
    """

    table_buffer = io.StringIO()
    with open(table_file, "r") as table:
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


def normalize_output_table(table_file, query_name, query_length, score_threshold, evalue_threshold):

    # note delimiter here: prehistoric style, several ws chars to seperate columns ...
    # additionally, will fail if hash used as separator --- see function
    # load_table_into_buffer
    try:
        df = pd.read_csv(
            table_file,
            sep="\s+",
            header=None,
            skip_blank_lines=True,
            comment="#",
            names=HMMER_TABLE_NAMES,
            usecols=HMMER_TABLE_USE_COLS,
            dtype=HMMER_TABLE_DTYPES
        )
    except ValueError:
        err_msg = (
            f"\nError parsing file: {table_file}\n"
            "Standard parsing with pandas.read_csv(... comment='#' ...) failed.\n"
            "Assuming that '#' was used as part of an identifier in the file.\n"
            "Loading table into buffer and restart...\n\n"
        )
        sys.stderr.write(err_msg)
        table_buffer = load_table_into_buffer(table_file)
        df = pd.read_csv(
            table_buffer,
            sep="\s+",
            header=None,
            names=HMMER_TABLE_NAMES,
            usecols=HMMER_TABLE_USE_COLS,
            dtype=HMMER_TABLE_DTYPES
        )

    if df.empty:
        # rare case
        return df

    assert df["query"].nunique() == 1
    assert df["query"].unique()[0] == query_name, f"Table file does not match motif: {query_name}"

    df["bit_score"] = df["bit_score"].round(1)

    df["query_length"] = query_length
    df["query_hit_pct"] = (df["query_hit_end"] - (df["query_hit_start"] - 1)) / df["query_length"] * 100
    df["query_hit_pct"] = df["query_hit_pct"].round(2)

    if score_threshold != -1 or evalue_threshold != -1:
        df["hit_hiq"] = 0
        select_above_score_t = df["bit_score"] > score_threshold
        if evalue_threshold == -1:
            select_below_evalue_t = df["evalue"].astype(float) < sys.maxsize
        else:
            select_below_evalue_t = df["evalue"].astype(float) < evalue_threshold
        select_is_hiq = select_above_score_t & select_below_evalue_t
        df.loc[select_is_hiq, "hit_hiq"] = 1
    else:
        df["hit_hiq"] = -1

    # for hits on the revcomp strand, switch start and end to create valid BED output;
    # if the motif were to be extracted, the coordinates still match but the sequence
    # would have to be revcomp'ed to (exactly) match the motif sequence
    revcomp_select = df["target_strand"] == "-"
    forward_cols = ["target_hit_start", "target_hit_end", "target_env_start", "target_env_end"]
    reverse_cols = ["target_hit_end", "target_hit_start", "target_env_end", "target_env_start"]
    df.loc[revcomp_select, forward_cols] = df.loc[revcomp_select, reverse_cols].values
    assert (df["target_hit_start"] < df["target_hit_end"]).all()

    # NB: we keep the sort order generated by nhmmer/HMMER

    return df


def aggregate_hits_per_target(motif_hits):

    records = []
    for trg, hits in motif_hits.groupby("target"):
        record = col.OrderedDict({
            "target": trg,
            "query": hits["query"].iloc[0],
            "target_length": hits.at[hits.index[0], "target_length"],
            "num_hits_total": hits.shape[0],
            "num_hits_hiq": hits.loc[hits["hit_hiq"] == 1, :].shape[0],
            "num_hits_loq": hits.loc[hits["hit_hiq"] == 0, :].shape[0],
        })
        if record["num_hits_hiq"] == 0:
            # if no high-quality hits were found,
            # add default empty fields
            record["pct_hits_hiq"] = 0
            record["num_bp_hiq"] = 0
            record["pct_bp_hiq"] = 0
            record["mean_pct_len_hiq"] = 0
            record["median_pct_len_hiq"] = 0
            record["mean_score_hiq"] = 0
            record["median_score_hiq"] = 0
            records.append(record)
            continue
        record["pct_hits_hiq"] = round(record["num_hits_hiq"] / record["num_hits_total"] * 100, 2)
        hiq_subset = hits.loc[hits["hit_hiq"] > 0, :].copy()
        # -1 in following calculation
        record["num_bp_hiq"] = (hiq_subset["target_hit_end"] - (hiq_subset["target_hit_start"] - 1)).sum()
        record["pct_bp_hiq"] = round(record["num_bp_hiq"] / record["target_length"] * 100, 2)
        record["mean_pct_len_hiq"] = round(hiq_subset["query_hit_pct"].mean(), 2)
        record["median_pct_len_hiq"] = round(hiq_subset["query_hit_pct"].median(), 2)
        record["mean_score_hiq"] = round(hiq_subset["bit_score"].mean(), 0)
        record["median_score_hiq"] = round(hiq_subset["bit_score"].median(), 0)
        records.append(record)

    agg = pd.DataFrame.from_records(records).fillna(0, inplace=False)
    agg["num_bp_hiq"] = agg["num_bp_hiq"].astype(int)

    assert all(column in agg.columns for column in HMMER_AGG_COLUMNS)
    assert all(column in HMMER_AGG_COLUMNS for column in agg.columns)

    return agg


def reduce_to_bed_like_table(motif_hits):

    bed_like = motif_hits[HMMER_BED_COLUMNS].copy()
    bed_like.sort_values(["target", "target_hit_start", "target_hit_end"], inplace=True)
    # BED format is zero-based, nhmmer output is not
    bed_like["target_hit_start"] -= 1

    select_hiq = bed_like["hit_hiq"] == 1

    return bed_like, select_hiq


def dump_file_to_disk(dataframe, file_path, prefix_header=False):

    file_path.parent.mkdir(exist_ok=True, parents=True)
    if prefix_header:
        # use xopen here to support compressed output
        # w/o manually scanning file extension etc.
        with xopen.xopen(file_path, "w") as dump:
            _ = dump.write("#")
            dataframe.to_csv(dump, sep="\t", header=True, index=False)
    else:
        dataframe.to_csv(file_path, sep="\t", header=True, index=False)
    # if that operation succeeded, a non-empty output
    # file was created, hence, delete a potential
    # EMPTY flag file that was created during a previous run
    remove_flag_file(file_path)
    return


def create_empty_flag_file(file_path):

    flag_file = get_empty_flag_file_path(file_path)
    with open(flag_file, "w"):
        pass
    return


def get_empty_flag_file_path(file_path):

    file_ext = file_path.suffix  # this has leading dot
    empty_flag = f"{file_ext}.EMPTY"
    flag_file = file_path.with_suffix(empty_flag)
    return flag_file


def remove_flag_file(file_path):

    flag_file = get_empty_flag_file_path(file_path)
    flag_file.unlink(missing_ok=True)
    return


def main():

    args = parse_command_line()

    query_name, query_length = get_motif_info(args.motif_file)

    motif_hits = normalize_output_table(
        args.hmmer_tblout,
        query_name, query_length,
        args.score_threshold, args.evalue_threshold
    )

    if not motif_hits.empty:

        dump_file_to_disk(motif_hits, args.norm_table)

        aggregated_hits = aggregate_hits_per_target(motif_hits)
        dump_file_to_disk(aggregated_hits, args.aggregated)

        bed_all, select_hiq = reduce_to_bed_like_table(motif_hits)
        dump_file_to_disk(bed_all, args.bed_all, True)
        if not select_hiq.any() and args.force_empty:
            empty_bed = pd.DataFrame(columns=HMMER_BED_COLUMNS)
            dump_file_to_disk(empty_bed, args.bed_hiq, True)
            if args.empty_flag:
                create_empty_flag_file(args.bed_hiq)
        else:
            dump_file_to_disk(bed_all.loc[select_hiq, :], args.bed_hiq, True)

    elif args.force_empty:
        dump_file_to_disk(motif_hits, args.norm_table)

        empty_agg = pd.DataFrame(columns=HMMER_AGG_COLUMNS)
        dump_file_to_disk(empty_agg, args.aggregated)

        empty_bed = pd.DataFrame(columns=HMMER_BED_COLUMNS)
        dump_file_to_disk(empty_bed, args.bed_all, True)
        dump_file_to_disk(empty_bed, args.bed_hiq, True)

        if args.empty_flag:
            for out_file in [args.norm_table, args.aggregated, args.bed_all, args.bed_hiq]:
                create_empty_flag_file(out_file)
    else:
        pass  # script will terminate peacefully w/o output

    return 0


if __name__ == "__main__":
    main()
