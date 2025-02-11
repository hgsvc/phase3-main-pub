#!/usr/bin/env python3

import argparse as argp
import pathlib as pl
import sys

import pandas as pd


# same thresholds as for chrY paper
Y_SPECIFIC_MOTIFS = ['DYZ1_Yq', 'DYZ18_Yq', 'DYZ3-sec_Ycentro']
UNSPECIFIC_MOTIFS = ['DYZ2_Yq']
UNSPECIFIC_THRESHOLD = 300

# see below: made script parameter
#PRIMARY_ALN_THRESHOLD_PCT = 90


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--contig-ref-align", "-a",
        type=lambda x: pl.Path(x).resolve(strict=True),
        nargs="+",
        dest="alignments",
        help="Contig to reference alignment files."
    )

    parser.add_argument(
        "--drop-alignments",
        "-da",
        type=int,
        default=2,
        dest="drop_alignments",
        help="Drop alignments of type 1/-1/2 (primary/inverted/secondary). Default: 2"
    )

    parser.add_argument(
        "--min-aln-threshold", "-maln",
        type=int,
        default=50,
        dest="min_aln_t",
        help="Min. alignments threshold. Default: 50"
    )

    parser.add_argument(
        "--motif-hits", "-m",
        type=lambda x: pl.Path(x).resolve(strict=True),
        nargs="+",
        dest="motifs",
        default=[],
        help="Motif hit summary files."
    )

    parser.add_argument(
        "--select-chrom", "-s",
        type=str,
        choices=["chrY", "chrX"],
        dest="select_chrom"
    )

    parser.add_argument(
        "--output", "-o",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="output",
        help="Output TSV table with selected contig summary info."
    )

    args = parser.parse_args()

    return args


def summarize_alignments(aln_subset, query, query_length):

    agg_aln = aln_subset.groupby(["target_name", "align_orient"])["align_matching"].sum()
    # collect info as simple list
    agg_aln = [(query, query_length, target, aln_bp, orient, "aln") for (target, orient), aln_bp in agg_aln.items()]
    return agg_aln


def read_alignment_file(file_path, select_chrom, drop_aln):

    aln = pd.read_csv(file_path, sep="\t", header=0)
    aln.drop(
        [
            "cg_cigar", "cs_tag_diffs", "s1_chaining_score", "s2_chaining_score",
            "rl_len_rep_seed_query", "cm_num_chain_miniz", "nn_ambig_bases"
        ],
        axis=1, inplace=True
    )
    if drop_aln != 0:
    # if set: drop alignments, i.e. most commonly secondary
        aln = aln.loc[aln["tp_align_type"] != drop_aln, :].copy()

    # aggregate over query name, skip all w/o any alignment to select_chrom
    query_infos = []
    for query, aligns in aln.groupby("query_name"):
        if select_chrom not in aligns["target_name"].values:
            continue
        query_length = aligns["query_length"].unique()[0]
        agg_aln = summarize_alignments(aligns, query, query_length)
        query_infos.extend(agg_aln)

    return query_infos


def read_motif_hit_file(file_path):

    hits = pd.read_csv(file_path, sep="\t", header=0)
    if hits.empty:
        return None
    if hits["query"].unique()[0] in Y_SPECIFIC_MOTIFS:
        # any high-quality hit counts
        subset = hits.loc[hits["num_hits_hiq"] > 0, :].copy()
        if subset.empty:
            return None
        query_infos = [(row.target, row.target_length, row.query, row.num_hits_hiq, 1, "hit") for row in subset.itertuples()]
    elif hits["query"].unique()[0] in UNSPECIFIC_MOTIFS:
        subset = hits.loc[hits["num_hits_hiq"] > UNSPECIFIC_THRESHOLD, :].copy()
        if subset.empty:
            return None
        query_infos = [(row.target, row.target_length, row.query, row.num_hits_hiq, -1, "hit") for row in subset.itertuples()]
    else:
        sys.stderr.write(f"\nWARNING: motif unknown - skipping ... {hits['query'].unique()[0]}\n")
        return None
    return query_infos


def check_alignment_ratio(contig_aln, select_chrom, min_aln_t):
    """For all contig that have not enough specific
    motif hits, check if the alignment ratio is above
    'min_aln_t' for the chromosome of interest.

    Args:
        contig_aln (_type_): _description_
    """

    aln_length = contig_aln.groupby(["tig_name", "ref_seq"])["num_bp_or_hits"].sum()
    tig_lengths = dict((row.tig_name, row.tig_length) for row in contig_aln.itertuples())

    selected = set()
    skipped = []
    for (tig, ref), tig_aln in aln_length.items():
        if ref != select_chrom:
            continue
        aln_ratio = round(tig_aln / tig_lengths[tig] * 100, 1)
        if aln_ratio > min_aln_t:
            selected.add(tig)
        else:
            skipped.append(
                f"# Skipped tig {tig} / {ref} - pct. aligned length: ~{aln_ratio}"
            )

    return selected, skipped


def main():

    args = parse_command_line()

    contig_infos = []
    for aln_file in args.alignments:
        alignments = read_alignment_file(aln_file, args.select_chrom, args.drop_alignments)
        contig_infos.extend(alignments)

    if args.select_chrom == "chrY":
        for motif_file in args.motifs:
            motif_hits = read_motif_hit_file(motif_file)
            if motif_hits is None:
                continue
            contig_infos.extend(motif_hits)

    df = pd.DataFrame.from_records(
        contig_infos,
        columns=["tig_name", "tig_length", "ref_seq", "num_bp_or_hits", "orient_or_specific", "info_type"]
    )

    # all contigs with specific or many unspecific motif hits will be selected by default
    select_has_hits = df["info_type"] == "hit"
    tigs_with_hits = set(df.loc[select_has_hits, "tig_name"].values)
    hit_subset = df.loc[df["tig_name"].isin(tigs_with_hits), :].copy()

    remaining_subset = df.loc[~df["tig_name"].isin(tigs_with_hits), :].copy()
    remaining_subset = remaining_subset.loc[remaining_subset["ref_seq"] == args.select_chrom, :].copy()

    high_aln_ratio, skipped_aln_ratio = check_alignment_ratio(remaining_subset, args.select_chrom, args.min_aln_t)
    select_high_aln = remaining_subset["tig_name"].isin(high_aln_ratio)
    aln_subset = remaining_subset.loc[select_high_aln, :].copy()

    selected_subset = pd.concat([hit_subset, aln_subset], axis=0, ignore_index=False)
    selected_subset.sort_values(["tig_length", "num_bp_or_hits"], ascending=False, inplace=True)

    args.output.parent.mkdir(exist_ok=True, parents=True)
    with open(args.output, "w") as table:
        table.write("\n".join(skipped_aln_ratio) + "\n")
        selected_subset.to_csv(table, sep="\t", header=True, index=False)


    return 0


if __name__ == "__main__":
    main()
