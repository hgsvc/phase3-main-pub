#!/usr/bin/env python3

import argparse as argp
import collections as col
import hashlib as hl
import pathlib as pl
import sys

import pandas as pd


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--contig-cov",
        "-c",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="contig_cov",
    )

    parser.add_argument(
        "--aln-unassign",
        "-u",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="aln_unassign",
    )

    parser.add_argument(
        "--aln-hap1",
        "-h1",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="aln_hap1",
    )

    parser.add_argument(
        "--aln-hap2",
        "-h2",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="aln_hap2",
    )

    parser.add_argument(
        "--aln-disconnect",
        "-d",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="aln_disconnect"
    )

    parser.add_argument(
        "--acrocentrics",
        "-a",
        type=str,
        nargs="+",
        dest="acrocentrics",
        default=["chr13", "chr14", "chr15", "chr21", "chr22"]
    )

    parser.add_argument(
        "--window-size",
        "-ws",
        type=int,
        dest="window_size",
        default=10000
    )

    parser.add_argument(
        "--out-target-view",
        "-ot",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="out_target"
    )

    parser.add_argument(
        "--out-query-view",
        "-oq",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="out_query"
    )

    parser.add_argument(
        "--out-stats",
        "-os",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="out_stats"
    )

    args = parser.parse_args()

    return args


def load_contig_coverages(norm_paf, acrocentrics):

    cov = pd.read_csv(
        norm_paf, sep="\t",
        header=[0, 1, 2],
        index_col=[0, 1, 2]
    )

    # drop chrM, chrX, chrY (for now)
    select_not_auto = cov.index.get_level_values("chrom").isin(["chrM", "chrX", "chrY"])
    cov = cov.loc[~select_not_auto, :].copy()

    acro_rows = cov.index.get_level_values("chrom").isin(acrocentrics)

    select_rdna = get_cov_selector("rdna")
    rdna_has_cov = (cov[select_rdna] > 0).any(axis=1)

    # create a boolean mask that indicates rows in the
    # coverage table that are likely uncovered by
    # hap1/hap2/unassign because they represent rDNA
    # arrays (separate in Verkko)
    rdna_mask = rdna_has_cov & acro_rows
    rdna_none = ~rdna_mask

    return cov, rdna_none


def aln_add_view_infos(row):
    target_view_aln = (
        f"QN:{row.query_name},"
        f"QL:{row.query_length},"
        f"OR:{row.align_orient},"
        f"QS:{row.query_start},"
        f"QE:{row.query_end},"
        f"MQ:{row.mapq}"
    )
    query_view_aln = (
        f"TN:{row.target_name},"
        f"QL:{row.query_length},"
        f"OR:{row.align_orient},"
        f"TS:{row.target_start},"
        f"TE:{row.target_end},"
        f"MQ:{row.mapq}"
    )
    return target_view_aln, query_view_aln


def load_alignments(hap1, hap2, unassign, disconnect):

    drop_columns = [
        "cg_cigar", "cs_tag_diffs", "s1_chaining_score", "s2_chaining_score",
        "rl_len_rep_seed_query", "cm_num_chain_miniz", "ms_dp_align_segment_max_score"
    ]

    aln_hap1 = pd.read_csv(hap1, sep="\t", header=0)
    aln_hap1.drop(drop_columns, axis=1, inplace=True)
    aln_hap1 = aln_hap1.loc[aln_hap1["tp_align_type"] == 1, :].copy()
    aln_hap1[["tview", "qview"]] = aln_hap1.apply(aln_add_view_infos, axis=1, result_type="expand")

    aln_hap2 = pd.read_csv(hap2, sep="\t", header=0)
    aln_hap2.drop(drop_columns, axis=1, inplace=True)
    aln_hap2 = aln_hap2.loc[aln_hap2["tp_align_type"] == 1, :].copy()
    aln_hap2[["tview", "qview"]] = aln_hap2.apply(aln_add_view_infos, axis=1, result_type="expand")

    aln_unassign = pd.read_csv(unassign, sep="\t", header=0)
    aln_unassign.drop(drop_columns, axis=1, inplace=True)
    aln_unassign = aln_unassign.loc[aln_unassign["tp_align_type"] == 1, :].copy()
    aln_unassign[["tview", "qview"]] = aln_unassign.apply(aln_add_view_infos, axis=1, result_type="expand")

    aln_disconnect = pd.read_csv(disconnect, sep="\t", header=0)
    aln_disconnect.drop(drop_columns, axis=1, inplace=True)
    aln_disconnect = aln_disconnect.loc[aln_disconnect["tp_align_type"] == 1, :].copy()
    aln_disconnect[["tview", "qview"]] = aln_disconnect.apply(aln_add_view_infos, axis=1, result_type="expand")

    alignments = {
        "hap1": aln_hap1,
        "hap2": aln_hap2,
        "unassigned": aln_unassign,
        "disconnected": aln_disconnect
    }
    return alignments


def get_cov_selector(asm_unit):
    selector = [
        (f"asm-{asm_unit}", "MQ00", "ctg_align_cov"),
        (f"asm-{asm_unit}", "MQ60", "ctg_align_cov"),
    ]
    return selector


def create_target_regions(tview_regions, tview_aln_records, issue_label, min_size=0):

    plain_regions = []
    for block_name, region_list in tview_regions.items():
        chrom = region_list[0][0]
        start = min(t[1] for t in region_list)
        end = max(t[2] for t in region_list)
        assert start < end
        if (end - start) < min_size:
            continue
        aln_records = tview_aln_records[block_name]
        plain_regions.append(
            (chrom, start, end, issue_label, block_name, aln_records)
        )
    plain_regions = pd.DataFrame.from_records(
        plain_regions, columns=["chrom", "start", "end", "name", "issue_id", "align_records"]
    )
    plain_regions.sort_values(["chrom", "start", "end"], inplace=True)

    return plain_regions


def create_query_regions(qview_nested_regions, qview_aln_records, issue_label, match_blocks):

    column_order = ["contig", "start", "end", "name", "issue_id", "align_records"]
    plain_regions = []
    for block_name, region_lists in qview_nested_regions.items():
        if block_name not in match_blocks:
            continue
        region_lists.reset_index(drop=False, inplace=True)
        region_lists["name"] = issue_label
        region_lists["issue_id"] = block_name
        region_lists["align_records"] = region_lists["issue_id"].replace(qview_aln_records)
        region_lists.rename(
            {
                "query_start": "start",
                "query_end": "end",
                "query_name": "contig"
            },
            inplace=True, axis=1
        )
        plain_regions.append(region_lists)
    plain_regions = pd.concat(plain_regions, axis=0, ignore_index=False)
    plain_regions = plain_regions[column_order]
    plain_regions = plain_regions.groupby(["contig", "start", "end", "name"]).agg(
        {
            "issue_id": lambda x: "@".join(x),
            "align_records": lambda x: "@".join(x)
        }
    )
    plain_regions.reset_index(drop=False, inplace=True)
    plain_regions.sort_values(["contig", "start", "end"], inplace=True)

    return plain_regions


def extract_diploid_regions(ctg_cov, rdna_free, alignments):

    select_h1 = get_cov_selector("hap1")
    select_h2 = get_cov_selector("hap2")
    select_un = get_cov_selector("unassigned")
    select_dis = get_cov_selector("disconnected")

    no_unassign = (ctg_cov[select_un] == 0).all(axis=1)
    no_disconn = (ctg_cov[select_dis] == 0).all(axis=1)
    any_h1 = (ctg_cov[select_h1] == 1).any(axis=1)
    any_h2 = (ctg_cov[select_h2] == 1).any(axis=1)

    selector = rdna_free & no_disconn & no_unassign & any_h1 & any_h2

    #diploid_regions = ctg_cov.loc[selector, :]
    #use_alignments = ["hap1", "hap2"]
    #label = "diploid"

    #tview_regions, qview_regions = process_selected_regions(
    #    diploid_regions, use_alignments, alignments, label
    #)

    return sum(selector)


def extract_dip_phasing_error(ctg_cov, rdna_free, alignments):

    select_h1 = get_cov_selector("hap1")
    select_h2 = get_cov_selector("hap2")
    select_un = get_cov_selector("unassigned")

    unassign = (ctg_cov[select_un] > 1).any(axis=1)
    no_h1 = (ctg_cov[select_h1] == 0).all(axis=1)
    no_h2 = (ctg_cov[select_h2] == 0).all(axis=1)

    selector = rdna_free & unassign & no_h1 & no_h2
    num_regions = sum(selector)

    if num_regions > 0:

        error_regions = ctg_cov.loc[selector, :]
        label = "dip_phasing_issue"
        use_alignments = ["unassigned"]

        tview_regions, qview_regions = process_selected_regions(
            error_regions, use_alignments, alignments, label
        )
    else:
        tview_regions = None
        qview_regions = None

    return tview_regions, qview_regions, num_regions


def extract_hap_phasing_error(ctg_cov, rdna_free, alignments, main, other):

    select_main = get_cov_selector(main)
    select_other = get_cov_selector(other)
    select_un = get_cov_selector("unassigned")

    unassign = (ctg_cov[select_un] > 0).any(axis=1)
    no_main = (ctg_cov[select_main] == 0).all(axis=1)
    some_other = (ctg_cov[select_other] > 0).any(axis=1)

    selector = rdna_free & unassign & no_main & some_other
    num_regions = sum(selector)

    if num_regions > 0:
        error_regions = ctg_cov.loc[selector, :]
        label = f"{main}_phasing_issue"
        use_alignments = ["unassigned"]

        tview_regions, qview_regions = process_selected_regions(
            error_regions, use_alignments, alignments, label
        )
    else:
        tview_regions = None
        qview_regions = None

    return tview_regions, qview_regions, num_regions


def extract_loh_regions(ctg_cov, rdna_free, alignments):

    select_h1 = get_cov_selector("hap1")
    select_h2 = get_cov_selector("hap2")
    select_un = get_cov_selector("unassigned")

    no_unassign = (ctg_cov[select_un] == 0).all(axis=1)
    hap_h1 = (ctg_cov[select_h1] == 1).all(axis=1)
    hap_h2 = (ctg_cov[select_h2] == 1).all(axis=1)

    selector = no_unassign & rdna_free & (hap_h1 ^ hap_h2)
    num_regions = sum(selector)

    if num_regions > 0:
        error_regions = ctg_cov.loc[selector, :]
        label = "LOH_misassm"
        use_alignments = ["hap1", "hap2"]

        tview_regions, qview_regions = process_selected_regions(
            error_regions, use_alignments, alignments, label, int(1e6)
        )
    else:
        tview_regions = None
        qview_regions = None

    return tview_regions, qview_regions, num_regions


def extract_shattered_regions(ctg_cov, rdna_free, alignments):

    select_h1 = get_cov_selector("hap1")
    select_h2 = get_cov_selector("hap2")
    select_un = get_cov_selector("unassigned")
    select_dis = get_cov_selector("disconnected")

    no_unassign = (ctg_cov[select_un] == 0).all(axis=1)
    disconn = (ctg_cov[select_dis] > 0).any(axis=1)
    no_h1 = (ctg_cov[select_h1] == 0).any(axis=1)
    no_h2 = (ctg_cov[select_h2] == 0).any(axis=1)

    selector = rdna_free & no_unassign & no_h1 & no_h2 & disconn
    num_regions = sum(selector)

    if num_regions > 0:
        error_regions = ctg_cov.loc[selector, :]
        label = "assm_broken"
        use_alignments = ["disconnected"]

        tview_regions, qview_regions = process_selected_regions(
            error_regions, use_alignments, alignments, label
        )
    else:
        tview_regions = None
        qview_regions = None

    return tview_regions, qview_regions, num_regions


def process_selected_regions(ctg_cov, use_alignments, alignments, label, min_size=0):

    annotated_regions = annotate_regions_with_alignments(ctg_cov, use_alignments, alignments)

    tview_regions = create_target_regions(
        annotated_regions["target_view_regions"],
        annotated_regions["target_view_aln_records"],
        label, min_size
    )

    qview_regions = create_query_regions(
        annotated_regions["query_view_regions"],
        annotated_regions["query_view_aln_records"],
        label, set(tview_regions["issue_id"].unique())
    )

    return tview_regions, qview_regions


def annotate_regions_with_alignments(issue_regions, use_alignments, alignments):

    # target view = BED-like track in reference coordinates
    target_view_regions = col.defaultdict(list)
    target_view_aln_records = dict()

    # query view = BED-like track in assembly coordinates
    query_view_regions = dict()
    query_view_aln_records = dict()

    for chrom, start, end in issue_regions.index:

        target_view_records = set()
        query_view_records = set()
        regions_by_query = dict()
        for aln_name in use_alignments:
            aln = alignments[aln_name]
            select_chrom = aln["target_name"] == chrom
            select_start = aln["target_end"] >= start
            select_end = aln["target_start"] < end

            sub = aln.loc[select_chrom & select_start & select_end, :]
            target_view_records = target_view_records.union(set(sub["tview"].values))
            query_view_records = query_view_records.union(set(sub["qview"].values))

            regions_by_query = sub.groupby("query_name")[["query_start", "query_end"]].agg({"query_start": min, "query_end": max})

        target_view_records = "|".join(sorted(target_view_records))
        query_view_records = "|".join(sorted(query_view_records))
        joined_name = target_view_records + "@" + query_view_records
        block_name = hl.md5(joined_name.encode("utf-8")).hexdigest()

        target_view_regions[block_name].append((chrom, start, end))
        target_view_aln_records[block_name] = target_view_records

        query_view_regions[block_name] = regions_by_query
        query_view_aln_records[block_name] = query_view_records

    annotated = {
        "target_view_regions": target_view_regions,
        "target_view_aln_records": target_view_aln_records,
        "query_view_regions": query_view_regions,
        "query_view_aln_records": query_view_aln_records
    }
    return annotated



def main():

    args = parse_command_line()

    cov, rdna_free = load_contig_coverages(args.contig_cov, args.acrocentrics)

    total_regions = cov.shape[0]
    other_regions = total_regions
    rdna_regions = (~rdna_free).sum()
    other_regions -= rdna_regions
    stats = [
        ("total_windows_autosomes", total_regions, int(total_regions * args.window_size)),
        ("rDNA_windows",rdna_regions, (rdna_regions * args.window_size))
    ]

    alignments = load_alignments(
        args.aln_hap1, args.aln_hap2,
        args.aln_unassign, args.aln_disconnect
    )

    target_view_out = []
    query_view_out = []

    num_dip_regions = extract_diploid_regions(cov, rdna_free, alignments)
    other_regions -= num_dip_regions
    stats.append(("diploid_regions", num_dip_regions, int(num_dip_regions * args.window_size)))

    tview_dip, qview_dip, num_dip = extract_dip_phasing_error(cov, rdna_free, alignments)
    other_regions -= num_dip
    if num_dip > 0:
        target_view_out.append(tview_dip)
        query_view_out.append(qview_dip)
    stats.append(("dip_phase_issue_windows", num_dip, int(num_dip * args.window_size)))

    tview_hap1, qview_hap1, num_hap1 = extract_hap_phasing_error(cov, rdna_free, alignments, "hap1", "hap2")
    other_regions -= num_hap1
    if num_hap1 > 0:
        target_view_out.append(tview_hap1)
        query_view_out.append(qview_hap1)
    stats.append(("hap1_phase_issue_windows", num_hap1, int(num_hap1 * args.window_size)))

    tview_hap2, qview_hap2, num_hap2 = extract_hap_phasing_error(cov, rdna_free, alignments, "hap2", "hap1")
    other_regions -= num_hap2
    if num_hap2 > 0:
        target_view_out.append(tview_hap2)
        query_view_out.append(qview_hap2)
    stats.append(("hap2_phase_issue_windows", num_hap2, int(num_hap2 * args.window_size)))

    tview_loh, qview_loh, num_loh = extract_loh_regions(cov, rdna_free, alignments)
    other_regions -= num_loh
    if num_loh > 0:
        target_view_out.append(tview_loh)
        query_view_out.append(qview_loh)
    stats.append(("loh_assm_issue_windows", num_loh, int(num_loh * args.window_size)))

    tview_broken, qview_broken, num_broken = extract_shattered_regions(cov, rdna_free, alignments)
    other_regions -= num_broken
    if num_broken > 0:
        target_view_out.append(tview_broken)
        query_view_out.append(qview_broken)
    stats.append(("misassm_broken_windows", num_broken, int(num_broken * args.window_size)))
    stats.append(("other_windows", other_regions, int(other_regions * args.window_size)))

    target_view_out = pd.concat(target_view_out, axis=0, ignore_index=False)
    target_view_out.sort_values(["chrom", "start"], inplace=True)
    args.out_target.parent.mkdir(exist_ok=True, parents=True)
    with open(args.out_target, "w") as tsv:
        _ = tsv.write("#")
        target_view_out.to_csv(tsv, sep="\t", header=True, index=False)

    query_view_out = pd.concat(query_view_out, axis=0, ignore_index=False)
    query_view_out.sort_values(["contig", "start"], inplace=True)
    args.out_query.parent.mkdir(exist_ok=True, parents=True)
    with open(args.out_query, "w") as tsv:
        _ = tsv.write("#")
        query_view_out.to_csv(tsv, sep="\t", header=True, index=False)

    stats = pd.DataFrame.from_records(stats, columns=["statistic", "num_windows", "num_bp"])
    stats["window_size"] = int(args.window_size)
    args.out_stats.parent.mkdir(exist_ok=True, parents=True)
    stats.to_csv(args.out_stats, sep="\t", header=True, index=False)

    return 0


if __name__ == "__main__":
    sys.exit(main())
