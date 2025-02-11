#!/usr/bin/env python3

import argparse as argp
import pathlib as pl
import collections as col
import itertools as itt

import pandas as pd
import xopen

# initial segments read from intersect file
Segment = col.namedtuple("Segment", ["seq", "start", "end", "label", "bpl_context", "crs_context", "src_idx", "dbg_idx"])

# final / merged regions
Region = col.namedtuple("Region", ["seq", "start", "end", "label", "size", "other_info", "bpl_blocks", "crs_blocks"])

# this is the same as in the script find_discontigs.py
GAP_LABELS = ["ALNGAP", "BRKSEQ", "UNALN"]
OVL_LABELS = ["BRKOVL"]
COMPLEX_LABELS = ["SCONTAIN", "XCONTAIN"]
ISSUE_LABELS = GAP_LABELS + OVL_LABELS + COMPLEX_LABELS


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--alignment-intersection", "-i",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="intersection",
        help="TSV table of intersection betwenn base-level and coarse-grained alignment."
    )

    parser.add_argument(
        "--known-ngaps", "-n",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="ngaps",
        help="BED-like file of known N gaps (unresolved) sequence in base-level coordinates."
    )

    parser.add_argument(
        "--ngap-label", "-l",
        type=str,
        default="NGAP",
        dest="ngap_label",
        help="Set the label for known N gaps. Default: NGAP"
    )

    parser.add_argument(
        "--output", "-o",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="output",
        help="Path to output table (TSV, BED-like)"
    )

    parser.add_argument(
        "--add-label-description", "-d",
        action="store_true",
        default=False,
        dest="add_description",
        help="Add label description to output table prefixed with '##'. Default: False"
    )

    args = parser.parse_args()

    return args


def initialize_aligned_segment(row, src_idx=None):

    global ISSUE_LABELS

    segment = None
    incomplete_aln_context = False
    if row.bpl_label in ISSUE_LABELS and row.crs_label in ISSUE_LABELS:
        # None segment triggers produce_flagged_segment()
        pass
    elif row.bpl_label in ISSUE_LABELS:
        # region aligned by coarse-grained aligner
        bpl_context = None
        crs_context = f"{row.crs_aln_context}::{row.crs_label}::{row.crs_blockid}"
        segment = Segment(row.crs_seq, row.crs_start, row.crs_end, "ALN", bpl_context, crs_context, src_idx, row.Index)
        incomplete_aln_context = True
    elif row.crs_label in ISSUE_LABELS:
        # region aligned by baselevel aligner
        bpl_context = f"{row.bpl_aln_context}::{row.bpl_label}::{row.bpl_blockid}"
        crs_context = None
        segment = Segment(row.bpl_seq, row.bpl_start, row.bpl_end, "ALN", bpl_context, crs_context, src_idx, row.Index)
        incomplete_aln_context = True
    elif row.bpl_label != row.crs_label:
        # in case of disagreemnt between the aligners w.r.t. to the target
        # sequence, simply build the segment that covers more sequence
        if row.bpl_size > row.crs_size:
            bpl_context = f"{row.bpl_aln_context}::{row.bpl_label}::{row.bpl_blockid}"
            crs_context = None
            segment = Segment(row.bpl_seq, row.bpl_start, row.bpl_end, "ALN", bpl_context, crs_context, src_idx, row.Index)
            incomplete_aln_context = True
        elif row.crs_size > row.bpl_size:
            bpl_context = None
            crs_context = f"{row.crs_aln_context}::{row.crs_label}::{row.crs_blockid}"
            segment = Segment(row.crs_seq, row.crs_start, row.crs_end, "ALN", bpl_context, crs_context, src_idx, row.Index)
            incomplete_aln_context = True
        else:
            #raise ValueError(f"Identical alignment size / different alignment targets: {row}")

            # empirically, that indeed does happen (rarely) but only with very short and
            # thus likely garbage assembled sequences. Since we only reach this code path
            # if the two aligners disagree on the correct placement, there is no other heuristic
            # to decide which alignment to favor. Hence, we assume that the base-level
            # aligner has stronger evidence for its placement of the sequence.
            bpl_context = f"{row.bpl_aln_context}::{row.bpl_label}::{row.bpl_blockid}"
            crs_context = None
            segment = Segment(row.bpl_seq, row.bpl_start, row.bpl_end, "ALN", bpl_context, crs_context, src_idx, row.Index)
            incomplete_aln_context = True
    else:
        # labels (= aligned target sequence) are identical
        assert row.bpl_seq == row.crs_seq
        seg_start = min(row.bpl_start, row.crs_start)
        seg_end = max(row.bpl_end, row.crs_end)
        bpl_context = f"{row.bpl_aln_context}::{row.bpl_label}::{row.bpl_blockid}"
        crs_context = f"{row.crs_aln_context}::{row.crs_label}::{row.crs_blockid}"
        segment = Segment(row.crs_seq, seg_start, seg_end, "ALN", bpl_context, crs_context, src_idx, row.Index)

    return segment, incomplete_aln_context


def produce_flagged_segment(row):
    """
    Produces a segment that indicates some form of issue
    """
    global GAP_LABELS
    global OVL_LABELS
    global ISSUE_LABELS
    flag_start = max(row.bpl_start, row.crs_start)
    flag_end = min(row.bpl_end, row.crs_end)

    # both are "proper" gaps
    both_gaps = row.bpl_label in GAP_LABELS and row.crs_label in GAP_LABELS
    # one is BRKOVL and the other indicates a gap, which is not a
    # gap in the alignment space, but definitely a discontinuity
    # in the assembly space.
    # If BRKOVL coincides with COMPLEX, error on the side of caution
    # and label the region as DISCON.
    is_discon = (
        row.bpl_label in GAP_LABELS and row.crs_label in OVL_LABELS
        or
        row.bpl_label in OVL_LABELS and row.crs_label in GAP_LABELS
        or
        row.bpl_label in OVL_LABELS and row.crs_label in OVL_LABELS
        or
        row.bpl_label in OVL_LABELS and row.crs_label in COMPLEX_LABELS
        or
        row.bpl_label in COMPLEX_LABELS and row.crs_label in OVL_LABELS
    )
    is_complex = (
        row.bpl_label in GAP_LABELS and row.crs_label in COMPLEX_LABELS
        or
        row.bpl_label in COMPLEX_LABELS and row.crs_label in GAP_LABELS
        or
        row.bpl_label in COMPLEX_LABELS and row.crs_label in COMPLEX_LABELS
    )

    if both_gaps:
        label = "ASMGAP"
    elif is_discon:
        label = "DISCON"
    elif is_complex:
        # logically, one or both must be complex
        label = "COMPLEX"
    else:
        raise RuntimeError(f"Invalid function call - not a flagged segment: {row}")
    bpl_context = f"{row.bpl_aln_context}::{row.bpl_label}::{row.bpl_blockid}"
    crs_context = f"{row.crs_aln_context}::{row.crs_label}::{row.crs_blockid}"
    segment = Segment(row.bpl_seq, flag_start, flag_end, label, bpl_context, crs_context, None, row.Index)
    return segment


def merge_context(context_infos):
    """This ugly function attempts to condense
    the alignment context info from both alignment types
    into one human-readable string. Candidate for breaking ...
    """

    spans = col.Counter()
    all_starts = col.defaultdict(list)
    all_ends = col.defaultdict(list)

    sample_info = set()
    block_ids = set()
    gap_contexts = set()
    abundance_counts = col.Counter()
    for context in context_infos:
        # key part NONE --- see function initialize_aligned_segment
        # the other: see N gap loading
        if context is None or context in ["KNOWN"]:
            continue
        try:
            src_context, other_seq, block_id = context.split("::")
        except ValueError:
            print(context)
            print(context_infos)
            raise
        block_ids.add(block_id)
        context_parts = src_context.split("|")
        sample_info.add(context_parts[0])
        other_seq = context_parts[1]
        if len(context_parts) < 4:
            # gaps are by definition w/o alignment information
            gap_contexts.add(src_context)
            abundance_counts[context_parts[0]] += 1
            abundance_counts[context_parts[0]] += 1
        else:
            other_start, other_end = context_parts[2].split("-")
            try:
                other_start = int(other_start)
                other_end = int(other_end)
            except ValueError:
                print(context_parts)
                raise
            assert other_start < other_end
            all_starts[other_seq].append(other_start)
            all_ends[other_seq].append(other_end)
            spans[other_seq] += (other_end - other_start)

    # select most likely correct other seq.
    try:
        select_seq = spans.most_common(1)[0][0]
    except IndexError:
        if len(gap_contexts) > 0:
            # record gap context
            select_seq = abundance_counts.most_common(1)[0][0]
            sample = abundance_counts.most_common(1)[0][0]
            start = 0
            end = 0
            assert len(block_ids) > 0
            all_blocks = "+".join(sorted(block_ids))
        else:
            select_seq = "UNK"
            sample = "UNK"
            all_blocks = "NO-BLOCK-INFO"
            start = 0
            end = 0
    else:
        start = min(all_starts[select_seq])
        end = max(all_ends[select_seq])
        all_blocks = "+".join(sorted(block_ids))
        assert len(sample_info) ==  1
        sample = sample_info.pop()
    return sample, select_seq, start, end, all_blocks


def merge_context_infos(bpl_contexts, crs_contexts):
    """Generate merged version of both alignment contexts
    and then decide which one to use (or combine) to
    roughly indicate the alignment context of the Region()
    (outside context)
    """

    bpl_sample, bpl_seq, bpl_start, bpl_end, bpl_blocks = merge_context(bpl_contexts)
    crs_sample, crs_seq, crs_start, crs_end, crs_blocks = merge_context(crs_contexts)

    if not bpl_sample == crs_sample and bpl_seq == crs_seq:
        raise ValueError(f"Incompatible contexts: {bpl_contexts} / {crs_contexts}")

    ignore_bpl = bpl_start == 0 and bpl_end == 0
    ignore_crs = crs_start == 0 and crs_end == 0

    if ignore_bpl and ignore_crs:
        info_field = f"{bpl_sample}|{bpl_seq}"
    elif ignore_bpl:
        info_field = f"{crs_sample}|{crs_seq}:{crs_start}-{crs_end}"
    elif ignore_crs:
        info_field = f"{bpl_sample}|{bpl_seq}:{bpl_start}-{bpl_end}"
    else:
        start = min(bpl_start, crs_start)
        end = max(bpl_end, crs_end)
        info_field = f"{bpl_sample}|{bpl_seq}:{start}-{end}"
    return info_field, bpl_blocks, crs_blocks


def get_seq(sequences):
    assert len(sequences) == 1
    return sequences.pop()


def get_label(ngap_label, labels):
    assert len(labels) > 0
    label = None
    if len(labels) > 1:
        assert ngap_label in labels
        label = ngap_label
    else:
        label = labels.pop()
    return label


def unify(segment_subset, ngap_label):
    """Unify is called on a subset of segments
    to perform the final coordinate and annotation
    merge. By construction, all segments must
    originate from the same sequence.
    """

    parts = []

    start = segment_subset[0].start
    end = segment_subset[0].end
    sequences = {segment_subset[0].seq}
    labels = {segment_subset[0].label}
    bpl_contexts = {segment_subset[0].bpl_context}
    crs_contexts = {segment_subset[0].crs_context}
    # add sentinel
    segment_subset.append(None)
    for a, b in itt.pairwise(segment_subset):
        try:
            if a.end >= b.start:
                end = b.end
                # in case in the next iteration,
                # there is no overlap, need to add
                # the b info already - consider case
                # of two overlapping segments
                bpl_contexts.add(b.bpl_context)
                crs_contexts.add(b.crs_context)
                labels.add(b.label)
                sequences.add(b.seq)
            elif a.end < b.start:
                bpl_contexts.add(a.bpl_context)
                crs_contexts.add(a.crs_context)
                sequences.add(a.seq)
                labels.add(a.label)

                info_field, bpl_blocks, crs_blocks = merge_context_infos(bpl_contexts, crs_contexts)
                parts.append(
                    Region(
                        get_seq(sequences), start, a.end, get_label(ngap_label, labels),
                        a.end - start, info_field, bpl_blocks, crs_blocks
                    )
                )
                bpl_contexts = {b.bpl_context}
                crs_contexts = {b.crs_context}
                labels = {b.label}
                sequences = {b.seq}
                start = b.start
                end = b.end
            else:
                raise
        except AttributeError:
            assert b is None
            break

    info_field, bpl_blocks, crs_blocks = merge_context_infos(bpl_contexts, crs_contexts)
    parts.append(
        Region(
            get_seq(sequences), start, end, get_label(ngap_label, labels),
            end-start, info_field, bpl_blocks, crs_blocks
        )
    )
    return parts


def build_segment_union(segments, ngap_label):
    """This function is not guaranteed
    to create a disjoined set of segments.
    This can be due to, e.g., enforced inclusion
    of N-gaps or overlapping alignments from different
    sequences.
    """

    total_union = []
    running_union = []
    # add sentinel to list of segments
    segments.append(None)
    # invariant in this loop:
    # only a is added to running_union
    for a, b in itt.pairwise(segments):
        try:
            same_seq = a.seq == b.seq
        except AttributeError:
            assert b is None
            # we reached the end of the list
            running_union.append(a)
            break
        if same_seq and a.label in ["ASMGAP", ngap_label] and b.label in ["ASMGAP", ngap_label]:
            running_union.append(a)
        elif same_seq and a.label == b.label:
            running_union.append(a)
        else:
            running_union.append(a)
            union_parts = unify(running_union, ngap_label)
            total_union.extend(union_parts)
            running_union = []
    if running_union:
        union_parts = unify(running_union, ngap_label)
        total_union.extend(union_parts)
    total_union = sorted(total_union, key=lambda r: (r.seq, r.start, r.end))

    return total_union


def read_ngap_segments(file_path, processed_sequences, ngap_label):

    ngaps = []
    idx = 0
    with xopen.xopen(file_path, "r") as table:
        for line in table:
            if line.startswith("#"):
                continue
            columns = line.strip().split()
            seq = columns[0]
            if seq not in processed_sequences:
                continue
            start = int(columns[1])
            end = int(columns[2])
            ngap = Segment(
                seq, start, end, ngap_label, "KNOWN", "KNOWN", None, f"ngap{idx}"
            )
            ngaps.append(ngap)
            idx += 1
    return ngaps


def get_intersection_header():
    """ Read intersect tables of the two alignment outputs:
    bpl = base-level alignments, e.g., produced with minimap2
    crs = coarse-grained alignments, e.g., produced with mashmap

    The prefixes bpl_ and crs_ are used throughout this script
    to distinguish these two components in the intersection tables.
    """

    commons = ["seq", "start", "end", "label", "aln_context", "size", "reflen", "blockid"]
    isect_columns = [f"bpl_{c}" for c in commons] + [f"crs_{c}" for c in commons] + ["overlap"]
    return isect_columns


def get_label_description_header(ngap_label, prefix=""):

    rows = []
    rows.append(
        f"{prefix}ALN: an aligned block supported by at least one aligner (base-level or coarse-grained)."
    )
    rows.append(
        (
            f"{prefix}ASMGAP: a gap, either in assembly-space (= broken sequence) or "
            "in alignment-space (= fragmented alignment / an alignment problem)."
        )
    )
    rows.append(
        (
            f"{prefix}DISCON: discontinuity in the assembly (two different sequences), "
            "but with an overlapping alignment in this region (= no simple gap)."
        )
    )
    rows.append(
        (
            f"{prefix}COMPLEX: alignments are fully contained in other alignments in this region. "
            "This can be an alignment artifact or indicate an assembly error and thus potentially "
            "a gap."
        )
    )
    rows.append(
        (
            f"{prefix}{ngap_label}: (if applicable) Known gaps due to unresolved sequence (N gaps). "
            "This annotation/labeling is only enforced for the respective coordinate space."
        )
    )

    return rows


def main():

    args = parse_command_line()

    NGAP_LABEL = args.ngap_label

    merged_segments = []
    active = False  # is there an active segment to be extended?
    current_segment = None
    incomplete_aln_context = False
    processed_sequences = set()  # used to limit ngaps to what is in the intersection

    df = pd.read_csv(args.intersection, sep="\t", header=None, names=get_intersection_header())

    for row in df.itertuples():
        processed_sequences.add(row.bpl_seq)
        # 0 --- check if we are inside an aligned block and
        # should just skip forward / skip this row.
        # This implements taking any alignment (base-level or coarse-grained)
        # as evidence that the assembly is contiguous and we thus
        # can skip over any gaps that the respective other aligner
        # is producing
        active = current_segment is not None  # just introduced as short-hand
        same_seq = active and row.bpl_seq == current_segment.seq
        bpl_inside = same_seq and row.bpl_end <= current_segment.end
        crs_inside = same_seq and row.crs_end <= current_segment.end
        skip_forward = bpl_inside and crs_inside
        if skip_forward:
            # in case of a disconnect between the aligners
            # in the start of the aligned block, record
            # a potentially missing context now
            contexts_available = row.bpl_label not in ISSUE_LABELS and row.crs_label not in ISSUE_LABELS
            # still skip if alignments to different sequences
            contexts_available &= (row.bpl_label == row.crs_label)
            if incomplete_aln_context and contexts_available:
                current_segment, incomplete_aln_context = initialize_aligned_segment(row, current_segment.dbg_idx)
            continue
        # 1 --- cannot skip, must process this row
        if active:
            # there is an active segment/alignment block,
            # so add this to the list of merged segments first
            merged_segments.append(current_segment)
            active = False
            incomplete_aln_context = False
            current_segment = None
        # initialize new current_segment with current row
        current_segment, incomplete_aln_context = initialize_aligned_segment(row)
        # 2 --- no alignment segment was initialized, this row
        # thus describes some form of gap/discontinuity
        if current_segment is None:
            # init flagged segment; note that active is also
            # still false
            flagged_segment = produce_flagged_segment(row)
            merged_segments.append(flagged_segment)
        # no else needed: init was succesful implies
        # that we are now in the next aligned block

    # iter complete, if there is still an active segment,
    # add it to the list
    if current_segment is not None:
        merged_segments.append(current_segment)

    # some post-processing; first, add known N-gaps
    ngaps = read_ngap_segments(args.ngaps, processed_sequences, NGAP_LABEL)
    merged_segments = sorted(merged_segments + ngaps, key=lambda seg: (seg.seq, seg.start, seg.end))

    # Next, build union of consecutive/overlapping segments.
    # This operation transforms the list from Segment() to Region()
    # and merges/simplifies the alignment contexts
    segment_union = build_segment_union(merged_segments, NGAP_LABEL)

    out = pd.DataFrame.from_records(
        segment_union, columns=["#seq", "start", "end", "label", "size", "align", "base_blocks", "coarse_blocks"]
    )
    args.output.parent.mkdir(exist_ok=True, parents=True)
    if args.add_description:
        desc_header = get_label_description_header(args.ngap_label, "## ")
        with xopen.xopen(args.output, "w") as bedlike:
            _ = bedlike.write("\n".join(desc_header) + "\n")
            out.to_csv(bedlike, sep="\t", header=True, index=False)
    else:
        out.to_csv(args.output, sep="\t", header=True, index=False)

    return 0


if __name__ == "__main__":
    main()
