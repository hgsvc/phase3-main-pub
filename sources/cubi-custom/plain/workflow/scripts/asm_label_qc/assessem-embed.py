#!/usr/bin/env python3

import argparse as argp
import collections as col
import pathlib as pl
import sys
import timeit

import pandas as pd
import pacmap
import numpy as np
import scipy.stats as stats


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--data-cache", "-dc",
        type=lambda x: pl.Path(x).resolve(strict=True),
        required=True,
        dest="data_cache",
        help="Path to existing data cache file (HDF format)."
    )

    parser.add_argument(
        "--bin-size", "-bs",
        type=int,
        default=int(1e4),
        dest="bin_size",
        help="Bin size in bp to aggregate cached data. Default: 10 000"
    )

    parser.add_argument(
        "--binned-data", "-ob", "--out-binned",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="binned_data",
        help="Path to output file holding the binned data."

    )

    parser.add_argument(
        "--transformed-data", "-ot", "--out-trans",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="transformed_data",
        help="Path to output file holding the transformed data."
    )

    parser.add_argument(
        "--skip-features", "-sf",
        type=str,
        nargs="*",
        dest="skip_features",
        default=None,
        help="Space-separated list of features to omit from embedding."
    )

    parser.add_argument(
        "--use-features", "-uf",
        type=str,
        nargs="*",
        dest="use_features",
        default=None,
        help="Space-separated list of features to use for embedding."
    )

    parser.add_argument(
        "--merge-features", "-mf",
        type=str,
        nargs="*",
        dest="merge_features",
        default=None,
        help="Space separated list of features to merge and average into single annotation"
    )

    parser.add_argument(
        "--verbose", "-vb",
        action="store_true",
        default=False,
        dest="verbose",
        help="Print logging and runtime information. Default: False"
    )

    args = parser.parse_args()


    return args


def print_runtime_summary(timer):

    for timed_op, timings in timer.items():
        delta = round(timings[3], 5)
        msg = (
            "\n==================\n"
            f"{timed_op}\n"
            f"{timings[0]} walltime: ~{delta} sec.\n"
            "==================\n"
        )
        sys.stdout.write(msg)
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


def build_final_data_track_list(track_labels, select_tracks, skip_tracks):
    """Validity of track labels has already been checked in the calling
    scope, so this function just does some set arithmetics and makes
    sure that no the user did not specify both "use" and "skip"
    """

    if select_tracks is None and skip_tracks is None:
        keep_tracks = sorted(track_labels)
    elif select_tracks is None:
        # skip tracks contains labels
        keep_tracks = sorted(track_labels - set(skip_tracks))
    elif skip_tracks is None:
        # select tracks contains labels
        pass
    else:
        # user wants to be funny
        info_msg = (
            "You specified both track labels to select and track "
            "labels to omit for the embedding process. That is "
            "illogical and only the 'select' labels will be honored "
            f"in the following: {select_tracks}"
        )
        sys.stderr.write(f"\n{info_msg}\n")
        keep_tracks = sorted(select_tracks)
    return keep_tracks



def load_cached_track_infos(data_cache, select_tracks, skip_tracks, merge_tracks):

    with pd.HDFStore(data_cache, "r") as hdf:
        try:
            track_infos = hdf["data_tracks"]
            track_labels = set(track_infos.index)
        except KeyError:
            err_msg = (
                f"The data cache file {data_cache} "
                "does not contain a track info item "
                "(key: 'data_tracks'). The cache seems "
                "to be invalid. Aborting."
            )
            raise RuntimeError(err_msg)

    unknown_select = set(select_tracks) - track_labels if select_tracks is not None else set()
    unknown_skip = set(skip_tracks) - track_labels if skip_tracks is not None else set()
    unknown_merge = set(merge_tracks) - track_labels if merge_tracks is not None else set()
    if not len(unknown_select) == len(unknown_skip) == len(unknown_merge) == 0:
        err_msg = (
            "At least one of the feature names / track labels you specified is not "
            "contained in the data cache file:\n"
            f"Unknown 'use-features': {unknown_select}\n"
            f"Unknown 'skip-features': {unknown_skip}\n"
            f"Unknown 'merge-features': {unknown_merge}\n"
        )
        raise ValueError(err_msg)

    # if set: determine which tracks to use for embedding
    keep_tracks = build_final_data_track_list(track_labels, select_tracks, skip_tracks)
    track_infos = track_infos.loc[keep_tracks, :].copy()

    if merge_tracks:
        # only binary or real-valued tracks can reasonbly be merged
        merge_types = set(track_infos.loc[merge_tracks, "score_type"].values)
        if len(merge_types) != 1 or merge_types.pop() not in ["binary", "quantitative"]:
            err_msg = (
                "Can only merge features/tracks of identical type, which must be "
                f"either 'binary' or 'quantitative'. Your selection: {merge_types}"
            )
            raise ValueError(err_msg)
        singletons = track_infos.index[~track_infos.index.isin(merge_tracks)]
        track_groups = [
            ([(row.Index, row.num_regions, row.score_type, row.data_range)], row.Index)
            for row in track_infos.loc[singletons, :].itertuples()
        ]
        track_groups.extend([
            (
                [
                    (row.Index, row.num_regions, None, None) for row in track_infos.loc[merge_tracks, :].itertuples()
                ],
                "merged")
        ])
    else:
        track_groups = [
            (
                [
                    (row.Index, row.num_regions, row.score_type, row.data_range)
                ],
                row.Index
            )
            for row in track_infos.itertuples()
        ]

    return track_groups


def load_region_track(data_cache, track_label, num_regions):

    with pd.HDFStore(data_cache, "r") as hdf:
        regions = hdf[f"tracks/{track_label}"]
    if not num_regions == regions.shape[0]:
        err_msg = (
            f"Size mismatch for number of regions: {track_label}\n"
            f"In metadata: {num_regions}\n"
            f"Loaded from cache: {regions.shape[0]}\n"
        )
        raise ValueError(err_msg)
    return regions


def load_cached_track_data(data_cache, track_infos, genome_size, bin_size, blunt_end):

    # TODO: track_info can be iterable
    # data_range defaults to float16 for merging tracks
    # num_regions can be collected when track info DF is loaded

    data_types = {"int8": np.int8, "uint8": np.uint8}
    if len(track_infos) > 1:
        # must be merge
        data_type = np.float16
    else:
        data_type = data_types[track_infos[0][3]]

    data = np.zeros(genome_size, dtype=data_type)
    last_score_type = None
    for (track_label, num_regions, score_type, _) in track_infos:
        regions = load_region_track(data_cache, track_label, num_regions)
        for row in regions.itertuples():
            data[row.start:row.end] += row.score
        data = data[:blunt_end]
        last_score_type = score_type

    if last_score_type is None:
        # implies merge data
        assert len(track_infos) > 1
        # compute average score
        data /= len(track_infos)

    data = data.reshape((-1, bin_size))
    if last_score_type is None:
        assert len(track_infos) > 1
        data = np.median(data, axis=1)
    elif last_score_type == "binary":
        data = (np.count_nonzero(data, axis=1)/bin_size).round(3).astype(np.float16)
    elif last_score_type == "categorical":
        data = np.apply_along_axis(lambda row: stats.mode(row).mode, 1, data)
    elif last_score_type == "quantitative":
        data = np.median(data, axis=1)
    else:
        raise ValueError(f"Unhandled score type: {last_score_type}")

    # the aggregate operations will result in floats
    # but high precision is pointless here
    data = data.astype(np.float16, casting="same_kind", copy=False)
    return data


def build_dataset(data_cache, bin_size, track_groups, timer):

    genome_size, _ = load_cached_genome_infos(data_cache)
    num_bins = genome_size // bin_size
    blunt_end = num_bins * bin_size

    num_tracks = len(track_groups)

    timer["alloc"] = [
        f"Allocating memory: {num_tracks}X{num_bins}",
        timeit.default_timer(),
        None, None
    ]
    dataset = np.zeros((num_tracks, num_bins), dtype=np.float16)
    timer["alloc"][2] = timeit.default_timer()
    timer["alloc"][3] = timer["alloc"][2] - timer["alloc"][1]

    label_order = []
    for row_num, (track_infos, track_label) in enumerate(track_groups, start=0):
        timer[f"load:{track_label}"] = [
            f"Loading track data: {[ti[0][0] for ti in track_infos]}",
            timeit.default_timer(),
            None, None
        ]
        dataset[row_num, :] = load_cached_track_data(data_cache, track_infos, genome_size, bin_size, blunt_end)
        label_order.append(track_label)
        timer[f"load:{track_label}"][2] = timeit.default_timer()
        timer[f"load:{track_label}"][3] = timer[f"load:{track_label}"][2] - timer[f"load:{track_label}"][1]

    dataset = dataset.transpose()
    dataset = pd.DataFrame(dataset, columns=label_order)
    return dataset


def compute_embedding_transformation(dataset):

    assert dataset.shape[0] > dataset.shape[1]

    # FOLLOWS default parameterization from quickstart
    # initializing the pacmap instance
    # Setting n_neighbors to "None" leads to a default choice shown below in "parameter" section
    embedding = pacmap.PaCMAP(n_components=2, n_neighbors=None, MN_ratio=0.5, FP_ratio=2.0)
    # fit the data (The index of transformed data corresponds to the index of the original data)
    transformed = embedding.fit_transform(dataset, init="pca")
    return transformed


def main():

    timer = col.OrderedDict()
    timer["assessem-embed"] = ["Total runtime", timeit.default_timer(), None, None]
    args = parse_command_line()

    # load cached track infos and limit to user selection, if applicable

    track_groups = load_cached_track_infos(
        args.data_cache, args.use_features,
        args.skip_features, args.merge_features
    )

    dataset = build_dataset(
        args.data_cache, args.bin_size,
        track_groups, timer
    )

    args.binned_data.parent.mkdir(exist_ok=True, parents=True)
    #dataset.to_feather(args.binned_data, compression='lz4')
    dataset.to_csv(args.binned_data, sep="\t", header=True, index=False)

    timer["embed"] = ["PaCMAP embedding", timeit.default_timer(), None, None]
    transformed = compute_embedding_transformation(dataset)
    timer["embed"][2] = timeit.default_timer()
    timer["embed"][3] = timer["embed"][2] - timer["embed"][1]

    args.transformed_data.parent.mkdir(exist_ok=True, parents=True)
    np.save(args.transformed_data, transformed, allow_pickle=False)

    timer["assessem-embed"][2] = timeit.default_timer()
    timer["assessem-embed"][3] = timer["assessem-embed"][2] - timer["assessem-embed"][1]

    if args.verbose:
        print_runtime_summary(timer)

    return 0


if __name__ == "__main__":
    main()
