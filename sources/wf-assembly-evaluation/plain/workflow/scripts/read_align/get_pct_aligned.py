#!/usr/bin/env python3

import argparse as argp
import pathlib as pl
import sys

import numpy as np
import pandas as pd
import pysam


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--input",
        "-i",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="input",
        help="Sorted and indexed BAM input file."
    )

    parser.add_argument(
        "--out-tsv",
        "-tsv",
        type=lambda x: pl.Path(x).resolve(),
        dest="out_tsv",
        help="Path to tab-separated output file listing all read alignment percentages."
    )

    parser.add_argument(
        "--out-hdf",
        "-hdf",
        type=lambda x: pl.Path(x).resolve(),
        dest="out_hdf",
        help="Path to HDF5 output file storing the actual coverage vectors."
    )

    parser.add_argument(
        "--min-pct-aligned",
        "-min-a",
        type=float,
        default=99.,
        dest="min_pct_aligned",
        help="Minimal percentage of read length that must be aligned to be counted. Default: 99%%"
    )

    parser.add_argument(
        "--min-read-length",
        "-min-l",
        type=int,
        default=1000,
        dest="min_read_length",
        help="Minimal (inferred) read length. Default: 1000"
    )

    parser.add_argument(
        "--io-threads",
        "-n",
        type=int,
        default=2,
        dest="io_threads",
        help="Number of threads for reading (decompressing) the BAM input file. Default: 2"
    )

    args = parser.parse_args()

    default_bai_file = args.input.with_suffix(".bam.bai")
    if not default_bai_file.is_file():
        sys.stderr.write(f"\nNo BAM index file found at path: {default_bai_file}\n")
        raise FileNotFoundError(default_bai_file)

    return args


def dump_read_table(table_file, align_records):

    file_mode = "w"
    add_header = True
    if table_file.is_file():
        file_mode = "a"
        add_header = False

    tsv_columns = ["seq_name", "read_name", "read_length", "aligned_pct", "mapq"]

    align_records = pd.DataFrame.from_records(
        align_records, columns=tsv_columns)

    align_records.to_csv(table_file, sep="\t", mode=file_mode, header=add_header, index=False)

    align_stats = align_records["aligned_pct"].describe()

    try:
        # these categories show up if the "object" datatype
        # is part of the statistics (happens for empty input)
        align_stats.drop(["top", "freq", "unique"], inplace=True)
    except KeyError:
        pass

    return align_stats


def dump_coverage(hdf_file, coverage, store_key):

    file_mode = "w"
    if hdf_file.is_file():
        file_mode = "a"

    coverage = pd.Series(coverage, dtype=int)
    with pd.HDFStore(hdf_file, mode=file_mode, complevel=9, complib="blosc") as hdf:
        hdf.put(store_key, coverage)

    cov_stats = coverage.describe(percentiles=[.25, .5, .75, .9, .95, .99])

    return cov_stats


def dump_statistics(hdf_file, align_stats):

    align_stats = pd.concat(align_stats, axis=1).transpose()
    align_stats.set_index("seq_name", inplace=True)
    align_stats.fillna(0., inplace=True)

    align_stats["seq_length"] = align_stats["seq_length"].astype(int)
    align_stats["alignments_total_count"] = align_stats["alignments_total_count"].astype(int)
    align_stats["reads_discarded_count"] = align_stats["reads_discarded_count"].astype(int)
    align_stats["alignments_selected_count"] = align_stats["alignments_selected_count"].astype(int)

    with pd.HDFStore(hdf_file, mode="a", complevel=9, complib="blosc") as hdf:
        hdf.put("statistics", align_stats)

    return


def main():

    args = parse_command_line()
    args.out_tsv.parent.mkdir(exist_ok=True, parents=True)
    args.out_hdf.parent.mkdir(exist_ok=True, parents=True)

    rename_align_stats = {
        "count": "alignments_total_count",
        "mean": "aligned_length_mean_pct",
        "std": "aligned_length_stddev_pct",
        "min": "aligned_length_min_pct",
        "25%": "aligned_length_q1_pct",
        "50%": "aligned_length_median_pct",
        "75%": "aligned_length_q3_pct",
        "max": "aligned_length_max_pct"
    }

    rename_cov_stats = {
        "count": "coverage_nonzero_count",
        "mean": "coverage_mean_abs",
        "std": "coverage_stddev_abs",
        "min": "coverage_min_abs",
        "25%": "coverage_q1_abs",
        "50%": "coverage_median_abs",
        "75%": "coverage_q3_abs",
        "90%": "coverage_90pct_abs",
        "95%": "coverage_95pct_abs",
        "99%": "coverage_99pct_abs",
        "max": "coverage_max_abs"
    }

    contig_align_stats = []

    # iterate over reference/target sequences from long to short
    with pysam.AlignmentFile(args.input, threads=args.io_threads) as bam:
        ref_seqs = sorted(
            [(int(length), name) for name, length in zip(bam.references, bam.lengths)],
            reverse=True
        )
        for ref_num, (ref_length, ref_name) in enumerate(ref_seqs, start=1):
            store_key = f"seq{ref_num}"
            align_stats = []
            coverage = np.zeros(ref_length, dtype=int)
            align_above_t = 0
            length_below_t = 0
            for record in bam.fetch(ref_name):
                if record.is_secondary or record.is_supplementary:
                    continue
                align_mapq = record.mapping_quality
                align_length = record.query_alignment_length
                read_length = record.infer_query_length()
                aligned_pct = round(align_length / read_length * 100, 1)
                align_stats.append((ref_name, record.query_name, read_length, aligned_pct, align_mapq))
                if read_length < args.min_read_length:
                    length_below_t += 1
                    continue
                if aligned_pct < args.min_pct_aligned:
                    continue
                coverage[record.reference_start:record.reference_end] += 1
                align_above_t += 1

            pct_stats = dump_read_table(args.out_tsv, align_stats)
            pct_stats.rename(rename_align_stats, inplace=True)
            pct_stats["store_key"] = store_key
            pct_stats["seq_name"] = ref_name
            pct_stats["seq_length"] = ref_length
            pct_stats["alignments_selected_count"] = align_above_t
            pct_stats["reads_discarded_count"] = length_below_t

            if align_above_t > 0:
                cov_stats = dump_coverage(args.out_hdf, coverage, store_key)
                cov_stats.rename(rename_cov_stats, inplace=True)
                cov_stats["coverage_nonzero_count"] = int((coverage > 0).sum())
                pct_stats = pd.concat([pct_stats, cov_stats], axis=0, ignore_index=False)

            contig_align_stats.append(pct_stats)

    dump_statistics(args.out_hdf, contig_align_stats)

    return 0


if __name__ == "__main__":
    main()
