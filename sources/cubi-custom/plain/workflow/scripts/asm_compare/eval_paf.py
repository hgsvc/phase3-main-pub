#!/usr/bin/env python3

import argparse as argp
import collections as col
import enum
import functools as fnt
import hashlib as hl
import math
import pathlib as pl
import re
import sys

import numpy as np
import numpy.ma as msk
import pandas as pd


class CigarOp(enum.Enum):
    IDENTITY = 0
    MISMATCH = 1
    MATCH = 2
    INSERTION = 3
    DELETION = 4
    HARDMASK = 5
    SOFTMASK = 6


class CigarNormalizer:

    def __init__(self):
        self.known_ops = {
            "=": CigarOp.IDENTITY,
            "X": CigarOp.MISMATCH,
            "M": CigarOp.MATCH,
            "I": CigarOp.INSERTION,
            "D": CigarOp.DELETION,
            "H": CigarOp.HARDMASK,
            "S": CigarOp.SOFTMASK
        }
        return

    def normalize(self, op_code):
        return self.known_ops[op_code]


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--paf-file", "-in",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="paf_file",
        help="PAF contig-to-contig alignment file (normalized)."
    )

    parser.add_argument(
        "--genome-size",
        "--fasta-index",
        "-s", "-f", "-gs",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="genome_size",
        default=None,
        help="Genome size file or FASTA index of the alignment target genome."
    )

    parser.add_argument(
        "--precision", "-prec",
        type=str,
        choices=["exact", "struct"],
        default="exact",
        dest="precision",
        help="Set desired precision: exact OR struct[ural]"
    )

    parser.add_argument(
        "--struct-err-geq", "-sst",
        type=int,
        default=50,
        dest="struct_err_geq",
        help=(
            "For precision 'struct[ural]', set minimum length to consider as error. "
            "In other words, deletions in the CIGAR string greater or equal "
            "(geq // >=) to this value will be recorded as errors. Default: 50"
        )
    )

    parser.add_argument(
        "--retain-mapq-geq", "-mq",
        "--keep-mapq-geq",
        type=int,
        default=1,
        dest="retain_mapq_geq",
        help=(
            "Retain only alignments that have a MAPQ greater or equal "
            "(geq // >=) to this value. Default: 1 (set to 0 to retain all)"
        )
    )

    parser.add_argument(
        "--retain-seqlen-geq", "-sl",
        "--keep-seqlen-geq",
        type=int,
        default=0,
        dest="retain_seqlen_geq",
        help=(
            "Retain only alignments where both the query and the target "
            "sequence have a length greater or equal (geq // >=) to this value. "
            "Default: 0"
        )
    )

    parser.add_argument(
        "--retain-alnlen-geq", "-aln",
        "--keep-alnlen-geq",
        type=int,
        default=10000,
        dest="retain_alnlen_geq",
        help=(
            "Retain only alignments where the length of the alignment is "
            "greater or equal (geq // >=) to this value in both the query "
            "and the target. Default: 10000"
        )
    )

    parser.add_argument(
        "--out-seq-stats", "-oss",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="out_seq_stats",
        default=None
    )

    parser.add_argument(
        "--out-aligned-regions", "-oar",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="out_aln_regions",
        default=None,
        help="Output TSV table with aligned regions including local QV estimate."
    )

    parser.add_argument(
        "--out-support-regions", "-osr",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="out_sup_regions",
        default=None,
        help=(
            "Output TSV table with support regions, i.e. regions identical "
            "under current parameterization between target and query."
        )
    )

    parser.add_argument(
        "--out-merged-regions", "-omr",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="out_mrg_regions",
        default=None,
        help=(
            "Output TSV table with merged regions, i.e. regions identical "
            "under current parameterization between target and query after "
            "merging overlapping regions (overlapping alignments)."
        )
    )

    parser.add_argument(
        "--split-seq-tag", "-tag",
        type=str,
        default=None,
        dest="split_seq_tag",
        help=(
            "Assume sequence names are tagged. Split off the tag using this "
            "character as separator. Default: None"
        )
    )


    parser.add_argument(
        "--verbose", "-vb",
        action="store_true",
        default=False,
        dest="verbose",
        help="Write process log to stdout. Default: False"
    )

    args = parser.parse_args()

    return args


def compute_qv(num_errors, ref_size):
    """_summary_

    Args:
        num_errors (_type_): _description_
        ref_size (_type_): _description_

    Returns:
        _type_: _description_
    """
    p = num_errors / ref_size
    try:
        q = -10 * math.log10(p)
    except ValueError:
        return 99
    return int(round(q, 0))


def read_genome_sizes(file_path, size_lower_bound):

    genome_sizes = dict()
    discarded = set()
    discard_log = []
    with open(file_path, "r") as listing:
        for line in listing:
            columns = line.strip().split()
            seq_name, seq_length = columns[:2]
            seq_length = int(seq_length)
            if seq_length <= size_lower_bound:
                discard_log.append(
                    (
                        "read_genome_sizes\tdiscard\t"
                        f"{seq_name}\t{seq_length}<={size_lower_bound}"
                    )
                )
                discarded.add(seq_name)
            genome_sizes[seq_name] = int(seq_length)
    return genome_sizes, discarded, discard_log


def filter_alignments(alignments, mapq_bound, size_bound, aligned_bound):

    discard_log = []
    discard = alignments["mapq"] < mapq_bound
    if discard.any():
        discard_log.append(
            (
                "filter_alignments\tdiscard\t"
                f"{discard.sum()}-aln-rows\tMAPQ<{mapq_bound}"
            )
        )
        # not sure what to record here ...
        #mapq_subset = alignments.loc[discard, :].copy()
        alignments = alignments.loc[~discard, :].copy()
    discard = (alignments["query_length"] < size_bound) | (alignments["target_length"] < size_bound)
    if discard.any():
        discard_log.append(
            (
                "filter_alignments\tdiscard\t"
                f"{discard.sum()}-aln-rows\t[QRY|TRG]LEN<{size_bound}"
            )
        )
        size_subset = alignments.loc[discard, :].copy()
        queries_dropped = size_subset["query_name"].nunique()
        targets_dropped = size_subset["target_name"].nunique()
        discard_log.append(
            (
                "filter_alignments\tdiscard\t"
                f"{queries_dropped}-num-seq\tQRYLEN<{size_bound}"
            )
        )
        discard_log.append(
            (
                "filter_alignments\tdiscard\t"
                f"{targets_dropped}-num-seq\tTRGLEN<{size_bound}"
            )
        )
        alignments = alignments.loc[~discard, :].copy()
    discard_query = (alignments["query_end"] - alignments["query_start"]) < aligned_bound
    discard_target = (alignments["target_end"] - alignments["target_start"]) < aligned_bound
    discard = discard_query | discard_target
    if discard.any():
        # note that his discards only alignments, not entire sequences
        discard_log.append(
            (
                "filter_alignments\tdiscard\t"
                f"{discard.sum()}-aln-rows\tALNLEN<={aligned_bound}"
            )
        )
        # not sure what to record here ...
        #mapq_subset = alignments.loc[discard, :].copy()
        alignments = alignments.loc[~discard, :].copy()

    return alignments, discard_log


def get_cigar_op_sets(precision):

    if precision == "exact":
        good_ops = set([CigarOp.IDENTITY])
        bad_ops = set([CigarOp.MISMATCH, CigarOp.DELETION])
        static_ops = set([CigarOp.INSERTION])
    elif precision == "struct":
        good_ops = set([CigarOp.IDENTITY, CigarOp.MISMATCH, CigarOp.MATCH])
        bad_ops = set([CigarOp.DELETION])
        static_ops = set([CigarOp.INSERTION])
    else:
        raise ValueError(f"Unknown level for precision parameter: {precision}")
    return good_ops, bad_ops, static_ops


def split_cigar_string(normalizer, size_dist, cigar_re, precision, struct_len_t, cigar_ops):
    """Function to parse and normalize series of CIGAR operations.
    If precision is set to structural, DELETIONs will be recoded as IDENTITY
    if the operation's step size is smaller / less than 'struct_len_t'

    Args:
        normalizer (_type_): _description_
        size_dist (_type_): _description_
        cigar_re (_type_): _description_
        precision (_type_): _description_
        struct_len_t (_type_): _description_
        cigar_ops (_type_): _description_

    Returns:
        _type_: _description_
    """

    op_series = []
    for cigar_op in cigar_re.finditer(cigar_ops):
        step = int(cigar_op.group("NUM"))
        norm_op = normalizer.normalize(cigar_op.group("OP"))
        # why recode only DELETION?
        # this is the only operation that is recorded as error
        # in the precision mode 'struct', see function get_cigar_op_sets(precision)
        if precision == "struct" and norm_op == CigarOp.DELETION:
            if step < struct_len_t:
                norm_op = CigarOp.IDENTITY
        op_series.append((step, norm_op))
        size_dist[norm_op].append(step)
    return op_series


def walk_cigar_ops(good_ops, bad_ops, static_ops, aln_start, aln_end, op_series):

    start = aln_start
    end = aln_start
    last_op = None
    for step, op_type in op_series:
        if op_type in good_ops:
            end += step
            last_op = op_type
        elif op_type in bad_ops:
            yield start, end
            # this steps over the
            # unsupported region
            start = end + step
            end = start
            last_op = op_type
        elif op_type in static_ops:
            last_op = op_type
        else:
            raise ValueError(f"Unsupported CIGAR operation: {op_type}")
    assert last_op in good_ops
    assert end == aln_end
    yield start, end

    return


def process_alignment_file(alignments, precision, struct_err_threshold):

    CigNorm = CigarNormalizer()
    cigar_op_sizes = col.defaultdict(list)
    CIGAR_OP_RE = re.compile("(?P<NUM>[0-9]+)(?P<OP>[MX\=DIHS]{1})")

    supported_regions = col.defaultdict(list)
    block_qvs = col.defaultdict(list)

    good_ops, bad_ops, static_ops = get_cigar_op_sets(precision)

    split_cigar = fnt.partial(
        split_cigar_string,
        CigNorm, cigar_op_sizes, CIGAR_OP_RE, precision, struct_err_threshold
    )
    walk_cigar = fnt.partial(walk_cigar_ops, good_ops, bad_ops, static_ops)

    for row in alignments.itertuples():
        row_id = hl.md5(f"{row.Index}{row.target_name}{row.target_start}{row.cg_cigar}".encode("utf-8")).hexdigest()
        op_series = split_cigar(row.cg_cigar)
        regions = [(start, end, row_id) for start, end in walk_cigar(row.target_start, row.target_end, op_series)]
        supported_bp = sum(region[1] - region[0] for region in regions)
        block_bp = row.target_end - row.target_start
        block_qv = compute_qv(block_bp - supported_bp, block_bp)
        block_qvs[row.target_name].append((row.target_start, row.target_end, row_id, block_qv))
        supported_regions[row.target_name].extend(regions)

    return cigar_op_sizes, supported_regions, block_qvs


def compute_sequence_statistics(seq_sizes, discarded, supported_regions, split_seq_tag):
    """Compute basic sequence statistics (number of support blocks, QV estimate etc.).
    For overlapping alignments, the computed QV over-estimates the quality
    and the statistics are hence also computed after merging overlapping support blocks.

    Args:
        seq_sizes (_type_): _description_
        discarded (_type_): _description_
        supported_regions (_type_): _description_
        split_seq_tag (_type_): _description_

    Returns:
        _type_: _description_
    """

    data_columns = [
        "support_regions", "support_length", "support_pct", "support_qv",
        "merged_sup_regions", "merged_sup_length", "merged_sup_pct", "merged_sup_qv"
    ]
    empty_row = [0] * len(data_columns)

    if split_seq_tag is not None:
        out_header = ["seq_name", "seq_tag", "seq_size"] + data_columns
        discarded = set(seq_name.rsplit(split_seq_tag, 1)[0] for seq_name in discarded)
    else:
        out_header = ["seq_name", "seq_size"] + data_columns

    seq_stats = []
    genome_merged_regions = col.defaultdict(list)
    for seq_name, seq_size in seq_sizes.items():
        if split_seq_tag is not None:
            out_name, seq_tag = seq_name.rsplit(split_seq_tag, 1)
            suffix = [out_name, seq_tag, seq_size]
        else:
            out_name = seq_name
            seq_tag = None
            suffix = [out_name, seq_size]
        try:
            seq_regions = supported_regions[seq_name]
            positions = np.zeros(seq_size, dtype=bool)
            total_support = 0
            for (start, end, _) in seq_regions:
                positions[start:end] = 1
                total_support += end - start

            # compute unmerged statistics
            num_regions = len(seq_regions)
            support_pct = round(total_support/seq_size * 100, 2)
            total_error = seq_size - total_support
            support_qv = compute_qv(total_error, seq_size)
            support_columns = [num_regions, total_support, support_pct, support_qv]

            # compute merged statistics
            merged_support = positions.sum()
            merged_regions = msk.clump_masked(
                msk.masked_array(np.zeros(seq_size, dtype=bool), mask=positions)
            )
            genome_merged_regions[seq_name].extend(
                [
                    (region_slice.start, region_slice.stop)
                    for region_slice in merged_regions
                ]
            )
            num_merged_regions = len(merged_regions)
            merged_sup_pct = round(merged_support/seq_size * 100, 2)
            merged_error = seq_size - merged_support
            merged_qv = compute_qv(merged_error, seq_size)
            merged_columns = [num_merged_regions, merged_support, merged_sup_pct, merged_qv]

            seq_stats.append(
                tuple(suffix + support_columns + merged_columns)
            )
        except KeyError:
            seq_stats.append(tuple(suffix + empty_row))

    seq_stats = pd.DataFrame.from_records(
        seq_stats, columns=out_header
    )
    seq_stats["discard_size"] = 0
    seq_stats.loc[seq_stats["seq_name"].isin(discarded), "discard_size"] = 1

    seq_stats.sort_values(["seq_name"], inplace=True)

    return seq_stats, genome_merged_regions


def create_region_table(regions, split_seq_tag, region_type):

    if region_type == "aligned":
        columns = ["start", "end", "region_id", "qv"]
    elif region_type == "merged":
        columns = ["start", "end"]
    elif region_type == "support":
        columns = ["start", "end", "region_id"]
    else:
        raise

    concat = []
    for seq_name, sub_regions in regions.items():
        df = pd.DataFrame.from_records(
            sub_regions, columns=columns
        )
        if split_seq_tag is not None:
            out_name = seq_name.rsplit(split_seq_tag, 1)[0]
        else:
            out_name = seq_name
        df["seq_name"] = out_name
        df = df[["seq_name"] + columns]
        concat.append(df)
    concat = pd.concat(concat, axis=0, ignore_index=False)
    concat.sort_values(["seq_name", "start", "end"], inplace=True)
    return concat


def log_process(log_lines, verbose):

    if log_lines and verbose:
        msg = "\n".join(log_lines) + "\n"
        sys.stdout.write(msg)
    return


def main():

    args = parse_command_line()

    gsize, discarded, log_lines = read_genome_sizes(args.genome_size, args.retain_seqlen_geq)
    log_process(log_lines, args.verbose)

    paf_aln = pd.read_csv(args.paf_file, sep="\t", header=0)
    paf_aln = paf_aln.loc[paf_aln["tp_align_type"] != 2, :].copy()

    paf_aln, log_lines = filter_alignments(
        paf_aln,
        args.retain_mapq_geq, args.retain_seqlen_geq, args.retain_alnlen_geq
    )
    log_process(log_lines, args.verbose)

    cgop_size_dist, supported_regions, aln_regions = process_alignment_file(paf_aln, args.precision, args.struct_err_geq)
    seq_stats, merged_regions = compute_sequence_statistics(gsize, discarded, supported_regions, args.split_seq_tag)

    if args.out_seq_stats is not None:
        args.out_seq_stats.parent.mkdir(parents=True, exist_ok=True)
        seq_stats.to_csv(args.out_seq_stats, sep="\t", header=True, index=False)

    if args.out_mrg_regions is not None:
        args.out_mrg_regions.parent.mkdir(parents=True, exist_ok=True)
        merged_regions = create_region_table(merged_regions, args.split_seq_tag, "merged")
        merged_regions.to_csv(args.out_mrg_regions, sep="\t", header=True, index=False)

    if args.out_aln_regions is not None:
        args.out_aln_regions.parent.mkdir(parents=True, exist_ok=True)
        aln_regions = create_region_table(aln_regions, args.split_seq_tag, "aligned")
        aln_regions.to_csv(args.out_aln_regions, sep="\t", header=True, index=False)

    if args.out_sup_regions is not None:
        args.out_sup_regions.parent.mkdir(parents=True, exist_ok=True)
        supported_regions = create_region_table(supported_regions, args.split_seq_tag, "support")
        supported_regions.to_csv(args.out_sup_regions, sep="\t", header=True, index=False)

    return 0


if __name__ == "__main__":
    main()
