#!/usr/bin/env python3

import argparse as argp
import collections as col
import pathlib as pl
import sys
import timeit

import pandas as pd
import numpy as np
import xopen

# these will be set as global
# constants depending on the command
# line argument "--allow-negative-scores"
_INT_SIZE = None
_CLIP_RANGE = None
_INT_REPR = None


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--genome-size", "-gs",
        type=lambda x: pl.Path(x).resolve(strict=True),
        default=None,
        dest="genome_size",
        help=(
            "2-column genome/contig size file: <NAME><TAB><SIZE>. "
            "Is required if the data cache file does not yet exist."
        )
    )

    parser.add_argument(
        "--data-cache", "-dc",
        type=lambda x: pl.Path(x).resolve(strict=False),
        required=True,
        dest="data_cache",
        help=(
            "Path to data cache file (HDF format). "
            "Will be created if it does not exist."
        )
    )

    parser.add_argument(
        "--allow-negative-scores", "-neg",
        action="store_true",
        default=False,
        dest="allow_negative",
        help=(
            "Allow negative scores in the data tracks in the range "
            "of -128 to 127. Default: False (= accept score within the "
            "range of 0 to 255)"
        )
    )

    parser.add_argument(
        "--overwrite-cached-entries", "--force",
        action="store_true",
        default=False,
        dest="force",
        help="Overwrite previously cached data tracks."
    )

    parser.add_argument(
        "--track-files", "-tf",
        type=lambda x: pl.Path(x).resolve(strict=True),
        nargs="+",
        dest="track_files",
        help="Path to data track files (TSV tables w/ header)."
    )

    parser.add_argument(
        "--track-labels", "-tl",
        type=str,
        nargs="+",
        dest="track_labels",
        help="One label per track file. Must be a valid Python identifier."
    )

    parser.add_argument(
        "--score-columns", "-sc",
        type=str,
        nargs="+",
        dest="score_columns",
        help=(
            "Identify score column in track file by name or by position (0-indexed). "
            "The special value 'binary' sets a fixed score of 1 for each region in "
            "the track file."
        )
    )

    parser.add_argument(
        "--threshold-categorical", "-tcat",
        type=int,
        default=10,
        dest="threshold_categorical",
        help=(
            "Assume the region scores represent a categorical variable "
            "if the number of unique values is below this threshold. "
            "Default: 10"
        )
    )

    parser.add_argument(
        "--verbose", "-vb",
        action="store_true",
        default=False,
        dest="verbose",
        help="Print logging information and runtime report. Default: False"
    )

    args = parser.parse_args()

    if not len(args.track_files) == len(args.track_labels) == len(args.score_columns):
        raise ValueError("Track files, labels and score columns must all have same length.")

    if not len(args.track_labels) == len(set(args.track_labels)):
        raise ValueError("Track labels must be unique")

    global _INT_SIZE
    global _CLIP_RANGE
    global _INT_REPR
    if args.allow_negative:
        _INT_SIZE = np.int8
        _CLIP_RANGE = (-128, 127)
        _INT_REPR = "int8"
    else:
        _INT_SIZE = np.uint8
        _CLIP_RANGE = (0, 255)
        _INT_REPR = "uint8"

    return args


def read_genome_size_file(table_file):

    if table_file is None:
        raise RuntimeError("Init of data cache file requires genome size file input")
    df = pd.read_csv(table_file, sep="\t", header=None, names=["sequence", "size"])
    df.sort_values(["size", "sequence"], ascending=[False, True], inplace=True)
    augmentation = []
    last_start = 0
    total_length = 0
    for pos, row in enumerate(df.itertuples(), start=1):
        start_coord = last_start
        end_coord = start_coord + row.size
        augmentation.append(
            (f"seq{pos}", start_coord, end_coord)
        )
        last_start += row.size
        total_length += row.size
    augmentation = pd.DataFrame.from_records(
        augmentation, index=df.index, columns=["alias", "start", "end"]
    )
    df = pd.concat([df, augmentation], axis=1, ignore_index=False)
    total_row = pd.DataFrame(
        [["genome", total_length, "genome", 0, total_length]],
        columns=df.columns
    )
    df = pd.concat([total_row, df], axis=0, ignore_index=False)
    df.reset_index(drop=True, inplace=True)
    return df


def init_data_cache_file(data_cache, genome_size, force):

    genome_size = read_genome_size_file(genome_size)
    data_cache.parent.mkdir(exist_ok=True, parents=True)

    with pd.HDFStore(data_cache, "w", complevel=1, complib="blosc") as hdf:
        hdf.put("genome_size", genome_size, format="fixed")
    return


def load_cached_genome_infos(data_cache):

    with pd.HDFStore(data_cache, "r") as hdf:
        try:
            genome_sizes = hdf["genome_size"]
        except KeyError:
            raise RuntimeError("Invalid data cache file - no entry 'genome_size'")
    total_length = genome_sizes.loc[genome_sizes["sequence"] == "genome", "size"].iloc[0]
    offsets = dict(
        (row.sequence, row.start) for row in genome_sizes.itertuples()
    )
    return total_length, offsets


def load_cached_track_infos(data_cache):

    with pd.HDFStore(data_cache, "r") as hdf:
        try:
            track_infos = hdf["data_tracks"]
            track_labels = set(track_infos.index)
        except KeyError:
            track_infos = None
            track_labels = set()
    return track_infos, track_labels


def parse_table_header(track_file, score_column):

    with xopen.xopen(track_file) as table:
        full_header = table.readline().strip().split()
    if score_column in full_header:
        use_columns = full_header[:3] + [score_column]
    else:
        try:
            col_idx = int(score_column)
            use_columns = full_header[:3] + full_header[col_idx]
        except ValueError:
            if score_column == "binary":
                use_columns = full_header[:3]
            else:
                raise ValueError(f"Cannot handle score column definition: {score_column}")
    return full_header, use_columns


def prepare_track_file_regions(track_file, score_column, offset_map, cat_threshold):

    full_header, load_columns = parse_table_header(track_file, score_column)
    regions = pd.read_csv(
        track_file, sep="\t", skiprows=1,
        header=None, names=full_header,
        usecols=load_columns
    )
    if len(load_columns) == 3:
        regions["score"] = 1
        new_names = ["seq", "start", "end"]
    else:
        new_names = ["seq", "start", "end", "score"]
    rename_columns = dict(
        (old_name, new_name) for old_name, new_name in zip(load_columns, new_names)
    )
    regions.rename(rename_columns, axis=1, inplace=True)

    # filter for regions that are in the offset_map
    # to ensure common set of regions across all annotations
    regions = regions.loc[regions["seq"].isin(offset_map), :].copy()

    regions["score"] = regions["score"].clip(
        lower=_CLIP_RANGE[0], upper=_CLIP_RANGE[1], inplace=False
    ).astype(_INT_SIZE)
    regions["coord_offset"] = regions["seq"].apply(lambda s: offset_map[s])
    regions["start"] += regions["coord_offset"]
    regions["end"] += regions["coord_offset"]
    regions = regions[["start", "end", "score"]].copy()

    uniq_scores = regions["score"].unique()
    if uniq_scores.min() == 0:
        # the annotations will be converted into
        # numpy arrays w/ a default of 0, hence
        # zero values should be ignored when determining
        # the score type
        uniq_scores = uniq_scores.size - 1
    else:
        uniq_scores = uniq_scores.size

    if uniq_scores < 2:
        score_type = "binary"
    elif uniq_scores < cat_threshold:
        score_type = "categorical"
    else:
        score_type = "quantitative"

    return regions, score_type


def add_track_to_data_cache(data_cache, track_file, track_label, regions, data_type, force):

    new_track = pd.DataFrame(
        [[data_type, regions.shape[0], _INT_REPR, str(track_file.name)]],
        index=[track_label],
        columns=["score_type", "num_regions", "data_range", "source_file"]
    )

    with pd.HDFStore(data_cache, "a", complevel=9, complib="blosc") as hdf:
        try:
            data_tracks = hdf["data_tracks"]
            if force:
                if track_label in data_tracks.index:
                    data_tracks.drop(track_label, axis=0, inplace=True)
            data_tracks = pd.concat([data_tracks, new_track], axis=0, ignore_index=False)
            data_tracks.sort_index(inplace=True)
            # no append due to data type size limitations
            hdf.put("data_tracks", data_tracks, format="fixed")
        except KeyError:
            hdf.put("data_tracks", new_track, format="fixed")
        hdf.put(f"tracks/{track_label}", regions, format="fixed")
    return


def print_runtime_summary(timer):

    for timed_op, timings in timer.items():
        delta = round(timings[3], 5)
        msg = (
            "\n==================\n"
            f"{timed_op}\n"
            f"{timings[0]} walltime: ~{delta} sec.\n"
            "====================\n"
        )
        sys.stdout.write(msg)
    return


def main():

    timer = col.OrderedDict()
    timer["assessem-cache"] = ["Total runtime", timeit.default_timer(), None, None]
    args = parse_command_line()
    if not args.data_cache.is_file():
        timer["init"] = [
            f"Init data cache file: {args.data_cache.name}",
            timeit.default_timer(),
            None, None
        ]
        init_data_cache_file(args.data_cache, args.genome_size, args.force)
        timer["init"][2] = timeit.default_timer()
        timer["init"][3] = timer["init"][2] - timer["init"][1]

    _, offset_map = load_cached_genome_infos(args.data_cache)
    _, cached_tracks = load_cached_track_infos(args.data_cache)

    zip_input = zip(args.track_files, args.track_labels, args.score_columns)
    for track_file, track_label, score_column in zip_input:
        if not track_label.isidentifier():
            if args.verbose:
                sys.stdout.write(f"\nTrack label is not a valid identifier: {track_label}")
                sys.stdout.write(f"\nSkipping over track file: {track_file}\n")
            continue
        if track_label in cached_tracks and not args.force:
            if args.verbose:
                sys.stdout.write(f"\nTrack label already in data cache: {track_label}")
                sys.stdout.write(f"\nSkipping over track file: {track_file}")
                sys.stdout.write("\n(set force overwrite to re-cache track)\n")
            continue

        timer[f"process:{track_label}"] = [
            f"Processing region track: {track_file.name}",
            timeit.default_timer(),
            None, None
        ]
        regions, data_type = prepare_track_file_regions(
            track_file, score_column, offset_map,
            args.threshold_categorical
        )
        timer[f"process:{track_label}"][2] = timeit.default_timer()
        timer[f"process:{track_label}"][3] = timer[f"process:{track_label}"][2] - timer[f"process:{track_label}"][1]

        timer[f"save:{track_label}"] = [
            f"Saving data track: {track_file.name}",
            timeit.default_timer(),
            None, None
        ]
        add_track_to_data_cache(args.data_cache, track_file, track_label, regions, data_type, args.force)
        timer[f"save:{track_label}"][2] = timeit.default_timer()
        timer[f"save:{track_label}"][3] = timer[f"save:{track_label}"][2] - timer[f"save:{track_label}"][1]

        cached_tracks.add(track_label)

    timer["assessem-cache"][2] = timeit.default_timer()
    timer["assessem-cache"][3] = timer["assessem-cache"][2] - timer["assessem-cache"][1]

    if args.verbose:
        print_runtime_summary(timer)

    return 0


if __name__ == "__main__":
    main()
