#!/usr/bin/env python3

import argparse as argp
import functools as fnt
import gzip
import pathlib as pl


import pandas as pd


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--regions", "-r",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="regions",
        help="Path to regions file in BED-like format (header required)."
    )

    parser.add_argument(
        "--read-cov", "-c",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="read_cov",
        help="Binned read coverage (HDF format)."
    )

    parser.add_argument(
        "--use-alignments", "-ua",
        type=int,
        choices=[0, 1],
        default=1,
        dest="aln_type",
        help=(
            "Annotate coverage using primary/1, supplementary/0 or secondary/2. "
            "Default: primary/1"
        )
    )

    parser.add_argument(
        "--use-mapq", "-uq",
        type=int,
        choices=[0, 60],
        default=0,
        dest="mapq",
        help=(
            "Annotate coverage using MAPQ 0 or MAPQ 60 alignments. "
            "Default: 0"
        )
    )

    parser.add_argument(
        "--table", "-t",
        "--output", "-o",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="table",
        help="Output table (TSV format)."
    )

    args = parser.parse_args()

    return args


def load_regions(regions_file):

    # TODO --- this script is currently executed in the snakemake
    # execution environment w/o the xopen module. That should be
    # added and the following simplified w/ xopen
    if regions_file.name.endswith("gz"):
        open_file = gzip.open
        open_mode = "rt"
    else:
        open_file = open
        open_mode = "r"

    with open_file(regions_file, open_mode) as bed_like:
        header = bed_like.readline().strip().strip("#").split()

    regions = pd.read_csv(
        regions_file, sep="\t", header=None,
        skiprows=[0], names=header
    )

    return regions


def get_contig_store_key(contig_map, contig_name):

    in_plain = contig_map.loc[contig_map["contig"] == contig_name, "store_key"]
    in_untagged = contig_map.loc[contig_map["untagged_contig"] == contig_name, "store_key"]
    if in_plain.empty and in_untagged.empty:
        raise ValueError(f"Cannot find store key for contig: {contig_name}")
    elif in_plain.empty:
        store_key = in_untagged.iloc[0]
    else:
        store_key = in_plain.iloc[0]
    return store_key


def annotate_regions_with_coverage(regions, cov_file, aln_type, mapq):
    """_summary_

    Args:
        regions (pandas.DataFrame): The regions to be annotated
        covdata (pathlib.Path): The HDF file containing the binned read coverage
        aln_type (int): Alignment type to select
        mapq (int): MAPQ value to select

    Raises:
        ValueError: _description_

    Returns:
        pandas.DataFrame: The annotated regions
    """

    with pd.HDFStore(cov_file, "r") as hdf:
        contig_map = hdf["contigs"]

        annotation = []
        ann_index = []
        idx = pd.IndexSlice
        last_contig = ""
        for region in regions.itertuples():
            # NB: the first column of the regions file could be
            # named chrom, contig, sequence etc., hence positional
            # access here to get the name: region[1] // idx 0 -> DataFrame.index
            contig_name = region[1]
            if contig_name != last_contig:
                contig_key = get_contig_store_key(contig_map, contig_name)
                covdata = hdf[contig_key]
                last_contig = contig_name
                covdata = covdata.xs((aln_type, mapq), axis=1, level=("aln_type", "mapq"))
                melted_labels = [f"{read}_{stat}" for read, stat in covdata.columns]
                covdata.columns = melted_labels
                _, start, end = covdata.index[0]
                win_size = end - start

            left = region.start // win_size * win_size
            # NB: index (= label) slicing in Pandas is right-inclusive,
            # hence one must *not* add 'win_size' to the right boundary
            # to have it included in the following slicing
            right = region.end // win_size * win_size
            windows = covdata.loc[idx[:, left:right, :], :]
            if windows.index.get_level_values("end").max() < region.end:
                # the last (incomplete) window is not included; this
                # just confirms that we are at the last region in that
                # case and we can ignore the missing window
                if region.Index != regions.index.max():
                    raise ValueError(f"Region out of bounds: {region} / {right}")
            windows = windows.median(axis=0).round(2)
            annotation.append(windows)
            ann_index.append(region.Index)

        annotation = pd.DataFrame(annotation, index=ann_index)
        regions = regions.merge(annotation, left_index=True, right_index=True, how="outer")

    return regions


def main():

    args = parse_command_line()

    regions = load_regions(args.regions)

    regions = annotate_regions_with_coverage(
        regions, args.read_cov,
        args.aln_type, args.mapq
    )

    args.table.parent.mkdir(exist_ok=True, parents=True)
    regions.to_csv(args.table, sep="\t", header=True, index=False)

    return 0


if __name__ == "__main__":
    main()
