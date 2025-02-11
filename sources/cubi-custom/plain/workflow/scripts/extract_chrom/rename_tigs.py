#!/usr/bin/env python3

import argparse as argp
import collections as col
import pathlib as pl

import dnaio
import pandas as pd


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--class-align", "-ca",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="cls_align",
    )

    parser.add_argument(
        "--seq-align", "-sa",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="seq_align",
    )

    parser.add_argument(
        "--input-class-order", "-io",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="class_order"
    )

    parser.add_argument(
        "--min-align-total",
        type=int,
        default=25000,
        dest="min_align_total"
    )

    parser.add_argument(
        "--fasta-in", "-fi",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="fasta_in"
    )

    parser.add_argument(
        "--fasta-out", "-fo",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="fasta_out",
        help="FASTA with renamed headers"
    )

    parser.add_argument(
        "--out-table", "-ot",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="out_table",
        help="Table with seq. class alignment summary."
    )

    parser.add_argument(
        "--out-class-align", "-oca",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="out_cls_align",
        help="Dump renamed alignment file (seq. classes to tigs)"
    )

    parser.add_argument(
        "--out-seq-align", "-osa",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="out_seq_align",
        help="Dump renamed alignment file (tigs to reference)"
    )

    args = parser.parse_args()

    return args


def read_alignment_table(file_path, min_align_size, keep_all=False, report_lost_contigs=None):
    """
    2024-05-06: added parameter report_lost_contigs for one edge case (CEPH/platinum pedegree)
    where one hifiasm assembly (sample: 200105-CEPH) has a Yq12/HET contig (~50 kbp, h1tg000857l.rev)
    that is end-to-end aligned in the sequence class to assembly alignment, but gap-aligned
    in the contig to sequence class case (for unclear reasons). Hence, this parameter is used
    to report contigs that are entirely lost after size thresholding in one of the two alignments
    and thus need to be removed in the other alignment table as well.
    """

    sample = file_path.name.split(".")[0]
    df = pd.read_csv(file_path, sep="\t", header=0)

    input_contigs = set()
    if report_lost_contigs is not None:
        input_contigs = set(df[report_lost_contigs].values)
    # this filtering step is the critical one; it virtually all cases
    # but the one described above in the doc string, this just gets rid
    # of spurious alignments
    df = df.loc[df["align_total"] > min_align_size, :].copy()
    retained_contigs = set()
    if report_lost_contigs is not None:
        retained_contigs = set(df[report_lost_contigs].values)

    lost_contigs = input_contigs - retained_contigs

    if not keep_all:
        # drop all "tag" headers
        drop_columns = [c for c in df.columns if len(c.split("_")[0]) == 2]
        df.drop(drop_columns, inplace=True, axis=1)

    return df, sample, lost_contigs


def get_seqclass_order(file_path):

    df = pd.read_csv(file_path, sep="\t", header=0)

    norm_names = {
        "other1": "OTH1",
        "HET1_centro": "HET1",
        "HET2_DYZ19": "HET2-DYZ19",
        "HET3_Yq": "HET3-Yq",
        "other2_DYZ18": "OTH2-DYZ18"
    }

    order_lut = dict()
    for order_num, row in enumerate(df.itertuples(), start=1):
        norm_name = norm_names.get(row.name, row.name)
        padded = f"{order_num:02}"
        order_lut[norm_name] = order_num, padded, f"{padded}|{norm_name}"

    df["name"] = df["name"].apply(lambda x: norm_names.get(x, x))
    df["padded"] = df["name"].apply(lambda x: order_lut[x][1])
    df["order_num"] = df["name"].apply(lambda x: order_lut[x][0])

    return order_lut, df


def pretty_print_length(length):

    if length > int(1e6):
        rounded_length = int(round(length, -6))
        if rounded_length > length:
            rounded_length -= int(1e6)
        rounded_length = rounded_length // int(1e6)
        prettified = f"{rounded_length}M"
    else:
        rounded_length = int(round(length, -3))
        while 1:
            if rounded_length < length:
                break
            rounded_length -= int(1e3)
        rounded_length = rounded_length // int(1e3)
        prettified = f"{rounded_length}k"

    assert rounded_length < length, f"{length} / {rounded_length}"

    return prettified


def annotate_sequences(sample, aligns, class_order):


    aligns["seq_class"] = aligns["query_name"].apply(lambda x: x.split("::")[0])
    aligns["order_num"] = aligns["seq_class"].apply(lambda x: class_order[x][0])

    group_tigs = col.defaultdict(list)
    tig_classes = dict()

    for tig, alignments in aligns.groupby("target_name"):
        ordered_aln = alignments.sort_values("order_num", inplace=False)
        start_class = ordered_aln["seq_class"].iloc[0]
        start_order = class_order[start_class][1]
        end_class = ordered_aln["seq_class"].iloc[-1]
        end_order = class_order[end_class][1]

        aligned_classes = set(ordered_aln["seq_class"].values)
        aligned_classes = sorted([class_order[cls][2] for cls in aligned_classes])

        tig_classes[tig] = aligned_classes

        tig_len = ordered_aln["target_length"].iloc[0]
        tig_len = pretty_print_length(tig_len)

        renamed = f"chrY|{start_order}.{end_order}|XX|{start_class}.{end_class}|1|{tig}|{tig_len}|{sample}"
        # if several tigs share an identical sequence class structure,
        # their relative order along the chromosome has to be determined
        # using the sequence/contig to reference alignment
        tig_rename_id = f"chrY|{start_order}.{end_order}|XX|{start_class}.{end_class}"
        group_tigs[tig_rename_id].append(
            (tig, renamed)
        )

    return group_tigs, tig_classes


def determine_tig_relative_order(seq_align, tig_infos):

    select_tigs = [t[0] for t in tig_infos]
    subset = seq_align.loc[seq_align["query_name"].isin(select_tigs), :].copy()

    aln_columns = ["query_length", "query_start", "query_end", "align_orient", "target_start", "target_end"]

    ordered_tigs = []
    for tig, alns in subset.groupby("query_name"):
        # find anchor - max align_matching
        aln_infos = alns.loc[alns["align_matching"].idxmax(), aln_columns]
        if aln_infos["align_orient"] > 0:
            remain_left = aln_infos["query_start"]
            remain_right = aln_infos["query_length"] - aln_infos["query_end"]
        else:
            remain_left = aln_infos["query_length"] - aln_infos["query_end"]
            remain_right = aln_infos["query_start"]

        # define an alignment bracket relative to the anchor
        # alignment and then select lowest start coordinate
        # in target/ref to determine relative contig order
        aln_bracket_left = aln_infos["target_start"] - remain_left
        aln_bracket_right = aln_infos["target_end"] + remain_right

        select_starts = subset["target_end"] > aln_bracket_left
        select_ends = subset["target_start"] < aln_bracket_right
        select_tig = subset["query_name"] == tig
        selector = select_starts & select_ends & select_tig

        min_target_start = subset.loc[selector, "target_start"].min()
        max_target_end = subset.loc[selector, "target_end"].max()
        ordered_tigs.append((min_target_start, max_target_end, tig))

    tig_order = dict()
    for pos, (start, end, tig) in enumerate(sorted(ordered_tigs), start=1):
        tig_order[tig] = f"{pos:02}", start, end

    return tig_order


def check_sequences_no_seq_class(fasta_input, name_lookup, seq_align_file, class_order_infos, fasta_output):

    no_sqcls_hits = []
    sequences = []
    table_out = []
    with dnaio.open(fasta_input) as fasta:
        for record in fasta:
            try:
                _ = name_lookup[record.name]
            except KeyError:
                new_name = None
                no_sqcls_hits.append((record.name, pretty_print_length(len(record.sequence))))
            sequences.append((record.name, record.sequence))

    if no_sqcls_hits:
        # assembled sequences that were not hit during the sequence class
        # alignment are commonly short-ish, and potentially garbage. Check
        # if they can be localized using an unfiltered alignment.
        raw_aln, sample, _ = read_alignment_table(seq_align_file, 0)
        tig_order = determine_tig_relative_order(raw_aln, no_sqcls_hits)
        for tig, (rel_pos, start, end) in tig_order.items():
            select_left = start >= class_order_infos["start"]
            select_right = start < class_order_infos["end"]
            selector = select_left & select_right
            if not selector.any():
                raise ValueError(f"no start selectable: {sample} / {tig}")
            start_class, start_pad_order = class_order_infos.loc[selector, ["name", "padded"]].values[0]

            select_left = end > class_order_infos["start"]
            select_right = end <= class_order_infos["end"]
            selector = select_left & select_right
            if not selector.any():
                raise ValueError(f"no end selectable: {sample} / {tig}")
            end_class, end_pad_order = class_order_infos.loc[selector, ["name", "padded"]].values[0]

            for unmap_tig, pp_len in no_sqcls_hits:
                if unmap_tig != tig:
                    continue
                new_name = f"chrY|{start_pad_order}.{end_pad_order}|{rel_pos}|{start_class}.{end_class}|2|{tig}|{pp_len}|{sample}"
                assert tig not in name_lookup
                name_lookup[tig] = new_name
                # NB: no hit OR hits that are too short/noisy
                table_out.append((tig, new_name, "NO-SEQCLASS-HIT"))
                break

    for tig, pp_len in no_sqcls_hits:
        if tig in name_lookup:
            continue
        new_name = f"chrY|99.99|99|UNK.UNK|3|{tig}|{pp_len}|{sample}"
        table_out.append((tig, new_name, "UNMAPPED"))
        assert tig not in name_lookup
        name_lookup[tig] = new_name

    fasta_output.parent.mkdir(exist_ok=True, parents=True)
    with dnaio.FastaWriter(fasta_output) as fasta:
        for old_name, sequence in sequences:
            renamed = name_lookup[old_name]
            fasta.write(renamed, sequence)

    return table_out, name_lookup


def dump_renamed_alignment_file(input_file, rename_lookup, column, output_file):

    df, _, _ = read_alignment_table(input_file, 0, True)
    df[column] = df[column].replace(rename_lookup)

    output_file.parent.mkdir(exist_ok=True, parents=True)
    df.to_csv(output_file, sep="\t", header=True, index=False)
    return



def main():

    args = parse_command_line()

    cls_align, sample, lost_cls_contigs = read_alignment_table(args.cls_align, args.min_align_total, report_lost_contigs="target_name")
    seq_align, sample2, lost_seq_contigs = read_alignment_table(args.seq_align, args.min_align_total, report_lost_contigs="query_name")
    assert sample == sample2

    lost_contigs = lost_cls_contigs.union(lost_seq_contigs)

    cls_align = cls_align.loc[~cls_align["target_name"].isin(lost_contigs), :].copy()
    seq_align = seq_align.loc[~seq_align["query_name"].isin(lost_contigs), :].copy()

    class_order, class_order_infos = get_seqclass_order(args.class_order)

    group_tigs, tig_classes = annotate_sequences(sample, cls_align, class_order)

    fixed_pos = "01"
    table_out = []
    rename_lookup = dict()
    for _, all_tigs in group_tigs.items():
        if len(all_tigs) == 1:
            tig, renamed = all_tigs[0]
            renamed = renamed.replace("XX", fixed_pos)
            aln_classes = tig_classes[tig]
            table_out.append(
                (tig, renamed, ",".join(aln_classes))
            )
            rename_lookup[tig] = renamed
        else:
            tig_order = determine_tig_relative_order(seq_align, all_tigs)
            for tig, renamed in all_tigs:
                tig_pos = tig_order[tig][0]
                aln_classes = tig_classes[tig]
                fmt_name = renamed.replace("XX", tig_pos)
                table_out.append(
                    (tig, fmt_name, ",".join(aln_classes))
                )
                rename_lookup[tig] = fmt_name

    # normalize alignment files and write to new output
    # check if unmapped sequences exist
    other_seq, rename_lookup = check_sequences_no_seq_class(
        args.fasta_in, rename_lookup, args.seq_align,
        class_order_infos, args.fasta_out
    )
    table_out.extend(other_seq)

    table_out = pd.DataFrame.from_records(
        table_out, columns=["old_name", "new_name", "aligned_seq_classes"]
    )
    table_out.sort_values("new_name", inplace=True)

    dump_renamed_alignment_file(
        args.cls_align, rename_lookup, "target_name", args.out_cls_align
    )

    dump_renamed_alignment_file(
        args.seq_align, rename_lookup, "query_name", args.out_seq_align
    )

    args.out_table.parent.mkdir(exist_ok=True, parents=True)
    table_out.to_csv(args.out_table, sep="\t", header=True, index=False)

    return 0


if __name__ == "__main__":
    main()
