#!/usr/bin/env python3

import argparse as argp
import collections as col
import pathlib as pl
import functools
import sys

import pandas as pd


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--gap-file", "-g",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="gap_file"
    )

    parser.add_argument(
        "--block-intersect", "-b",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="block_intersect"
    )

    parser.add_argument(
        "--nucfreq-regions", "-n",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="nucfreq_regions"
    )

    parser.add_argument(
        "--flagger-regions", "-f",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="flagger_regions"
    )

    parser.add_argument(
        "--chrom-assign", "-c",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="chrom_assign"
    )

    parser.add_argument(
        "--karyotype", "-k",
        type=str,
        dest="karyotype"
    )

    parser.add_argument(
        "--sample", "-s",
        type=str,
        dest="sample"
    )

    parser.add_argument(
        "--assembly", "-a",
        type=str,
        dest="assembly"
    )

    parser.add_argument(
        "--out-table", "-ot",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="out_table",
        default=None
    )

    parser.add_argument(
        "--out-summary", "-os",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="out_summary",
        default=None
    )

    args = parser.parse_args()

    return args


def load_gap_ids(file_path, karyotype):

    with open(file_path) as listing:
        header = listing.readline().strip().split()
        assert header[0].startswith("#")
        header[0] = "chrom"
        assert header[3] == "name"
        header[3] = "gap_id"
        header[1] = "gap_start"
        header[2] = "gap_end"

    gaps = pd.read_csv(
        file_path, comment="#", sep="\t",
        header=None, names=header,
        usecols=header[:-2]
    )

    # sanitize inputs
    drop_rows = gaps["gap_start"] == gaps["gap_end"]
    if drop_rows.any():
        sys.stderr.write(f"\nWARNING:\nDropping invalid rows / gaps of size 0 (start == end): {drop_rows.sum()}\n\n")
        gaps = gaps.loc[~drop_rows, :].copy()
        gaps.reset_index(inplace=True, drop=True)

    gaps = gaps[["chrom", "gap_start", "gap_end", "gap_id"]].copy()
    if karyotype == "female":
        gaps = gaps.loc[gaps["chrom"] != "chrY", :].copy()
    if karyotype == "male":
        gaps = gaps.loc[gaps["chrom"] != "chrX", :].copy()

    return gaps


def load_chrom_assign(file_path):

    df = pd.read_csv(file_path, sep="\t", comment="#", header=0)

    chrom_assign = dict()

    for query, assignments in df.groupby("query_name"):
        max_match = assignments["align_matching"].idxmax()
        chrom = assignments.at[max_match, "target_name"]
        chrom_assign[query] = chrom

    return chrom_assign


def read_flagger_regions(file_path):

    df = pd.read_csv(file_path, sep="\t", header=0)
    df.rename({"#contig": "contig"}, axis=1, inplace=True)
    df.drop(["score", "asm_unit"], axis=1, inplace=True)
    return df


def read_nucfreq_regions(file_path):

    df = pd.read_csv(file_path, sep="\t", header=0)
    df.rename({"#contig": "contig"}, axis=1, inplace=True)
    df.drop(["hifi_median_cov", "hifi_pct_median_cov", "ont_median_cov", "ont_pct_median_cov"], axis=1, inplace=True)
    return df


def summarize_regions(regions, contig, start, end):

    select_ctg = regions["contig"] == contig
    select_start = regions["start"] < end
    select_end = regions["end"] > start

    selector = select_ctg & select_start & select_end
    if not selector.any():
        summary = [0, 0]
    else:
        sub = regions.loc[selector, :].copy()
        first = sub.index[0]
        last = sub.index[-1]
        # the following ensures that only the aligned
        # region is considered and not flagged regions
        # beyond its limits
        sub.loc[first, "start"] = max(start, sub.at[first, "start"])
        assert sub.at[first, "end"] > start
        sub.loc[last, "end"] = min(end, sub.at[last, "end"])
        assert sub.at[last, "start"] < end
        total = (sub["end"] - sub["start"]).sum()
        pct = round(total / (end - start) * 100, 5)
        summary = [total, pct]
    return summary


@functools.lru_cache(2048)
def split_aln_context(align_context):

    _, aln_context = align_context.split("|")
    asm_ctg, asm_coord = aln_context.split(":")
    asm_start, asm_end = asm_coord.split("-")
    asm_start = int(asm_start)
    asm_end = int(asm_end)
    return asm_ctg, asm_start, asm_end



def process_gap_alnblock_table(sample, asm_unit, file_path, flagger_file, nucfreq_file, gaps, chrom_assign):

    df = pd.read_csv(file_path, sep="\t", header=0)
    df.drop(["aln_base_block", "aln_coarse_block"], axis=1, inplace=True)

    flagger = read_flagger_regions(flagger_file)
    nucfreq = read_nucfreq_regions(nucfreq_file)

    gap_status = []

    for row in gaps.itertuples():
        gap_id = row.gap_id
        gap_length = row.gap_end - row.gap_start
        assert gap_length > 0, row
        gap_chrom = row.chrom
        sub = df.loc[df["gap_id"] == gap_id, :]
        if sub.empty or (sub["aln_label"] == "ASMGAP").all():
            status_info = {
                "gap_id": gap_id, "gap_length": gap_length,
                "sample": sample, "asm_unit": asm_unit,
                "seq": "unknown", "start": -1, "end": -1,
                "overlap_pct": 0,
                "align_status": "open",
                "flagger_bp": 0, "flagger_pct": 0,
                "nucfreq_bp": 0, "nucfreq_pct": 0
            }
            gap_status.append(status_info)
        elif (sub["aln_label"] == "ALN").all() and (sub["overlap_pct"] > 99.9).all():

            match_found = False
            for row in sub.itertuples():
                ctg, start, end = split_aln_context(row.aln_info)
                try:
                    ctg_chrom = chrom_assign[ctg]
                except KeyError:
                    # example case: NA19347.vrk-ps-sseq // haplotype1-0000060
                    # a ~16 kbp sequence that only aligns with MAPQ 0 and is
                    # hence disregarded for the chromosome assignment
                    # --- assume this to be the standard case and skip this
                    continue
                if ctg_chrom != gap_chrom:
                    continue
                flagger_summary = summarize_regions(flagger, ctg, start, end)
                nucfreq_summary = summarize_regions(nucfreq, ctg, start, end)
                status_info = {
                    "gap_id": gap_id, "gap_length": gap_length,
                    "sample": sample, "asm_unit": asm_unit,
                    "seq": ctg, "start": start, "end": end,
                    "overlap_pct": row.overlap_pct,
                    "align_status": "covered",
                    "flagger_bp": flagger_summary[0], "flagger_pct": flagger_summary[1],
                    "nucfreq_bp": nucfreq_summary[0], "nucfreq_pct": nucfreq_summary[1]
                }
                gap_status.append(status_info)
                match_found = True

            if not match_found:
                status_info = {
                    "gap_id": gap_id, "gap_length": gap_length,
                    "sample": sample, "asm_unit": asm_unit,
                    "seq": "unknown", "start": -1, "end": -1,
                    "overlap_pct": 0,
                    "align_status": "open",
                    "flagger_bp": 0, "flagger_pct": 0,
                    "nucfreq_bp": 0, "nucfreq_pct": 0
                }
                gap_status.append(status_info)

        else:
            if (sub["aln_label"] == "ALN").any():
                match_found = False
                for row in sub.itertuples():
                    if row.aln_label != "ALN":
                        continue
                    ctg, start, end = split_aln_context(row.aln_info)
                    try:
                        ctg_chrom = chrom_assign[ctg]
                    except KeyError:
                        # example case: NA19347.vrk-ps-sseq // haplotype1-0000060
                        # a ~16 kbp sequence that only aligns with MAPQ 0 and is
                        # hence disregarded for the chromosome assignment
                        # --- assume this to be the standard case and skip this
                        continue
                    if ctg_chrom != gap_chrom:
                        continue
                    flagger_summary = summarize_regions(flagger, ctg, start, end)
                    nucfreq_summary = summarize_regions(nucfreq, ctg, start, end)
                    status_label = "partial"
                    if row.overlap_pct > 99.9:
                        status_label = "covered"
                    status_info = {
                        "gap_id": gap_id, "gap_length": gap_length,
                        "sample": sample, "asm_unit": asm_unit,
                        "seq": ctg, "start": start, "end": end,
                        "overlap_pct": row.overlap_pct,
                        "align_status": status_label,
                        "flagger_bp": flagger_summary[0], "flagger_pct": flagger_summary[1],
                        "nucfreq_bp": nucfreq_summary[0], "nucfreq_pct": nucfreq_summary[1]
                    }
                    gap_status.append(status_info)
                    match_found = True

                if not match_found:
                    status_info = {
                        "gap_id": gap_id, "gap_length": gap_length,
                        "sample": sample, "asm_unit": asm_unit,
                        "seq": "unknown", "start": -1, "end": -1,
                        "overlap_pct": 0,
                        "align_status": "open",
                        "flagger_bp": 0, "flagger_pct": 0,
                        "nucfreq_bp": 0, "nucfreq_pct": 0
                    }
                    gap_status.append(status_info)

            else:
                status_info = {
                    "gap_id": gap_id, "gap_length": gap_length,
                    "sample": sample, "asm_unit": asm_unit,
                    "seq": "unknown", "start": -1, "end": -1,
                    "overlap_pct": sub["overlap_pct"].sum(),
                    "align_status": "complex",
                    "flagger_bp": 0, "flagger_pct": 0,
                    "nucfreq_bp": 0, "nucfreq_pct": 0
                }
                gap_status.append(status_info)

    gap_status = pd.DataFrame.from_records(gap_status)
    return gap_status


def summarize_status(gap_table):

    thresholds = [101, 10, 1, 0.1]
    labels = ["closed_at_err_any", "closed_at_err_lst_10pct", "closed_at_err_lst_1pct", "closed_at_err_lst_01pct"]
    assert len(thresholds) == len(labels)

    summary = []
    # TODO: these count stats are discarded
    stats = col.Counter()

    for gap_id, gap_status in gap_table.groupby("gap_id"):
        if gap_status.shape[0] > 1:
            if gap_status["align_status"].nunique() > 1:
                if (gap_status["align_status"] == "covered").any():
                    clean_status = gap_status.loc[gap_status["align_status"] == "covered", :].copy()
                    assert clean_status.shape[0] == 1
                else:
                    raise RuntimeError(gap_status)
            else:
                clean_status = gap_status.copy()
            seqs = ";".join(sorted(clean_status["seq"].values))
            aln_label = clean_status["align_status"].iloc[0]
            if aln_label == "covered":
                flagger_pct = clean_status["flagger_pct"].iloc[0]
                nucfreq_pct = clean_status["nucfreq_pct"].iloc[0]

                last_label = None
                for t, l in zip(thresholds, labels):
                    if flagger_pct < t and nucfreq_pct < t:
                        stats[l] += 1
                        last_label = l
            else:
                last_label = "not_closed"
            summary.append(
                (
                    gap_id, clean_status["sample"].iloc[0], clean_status["asm_unit"].iloc[0],
                    seqs, aln_label, last_label, clean_status["gap_length"].iloc[0]
                )
            )
            stats[aln_label] += 1
        else:
            if gap_status["align_status"].iloc[0] == "covered":
                flagger_pct = gap_status["flagger_pct"].iloc[0]
                nucfreq_pct = gap_status["nucfreq_pct"].iloc[0]

                last_label = None
                for t, l in zip(thresholds, labels):
                    if flagger_pct < t and nucfreq_pct < t:
                        stats[l] += 1
                        last_label = l

                summary.append(
                    (
                        gap_id, gap_status["sample"].iloc[0], gap_status["asm_unit"].iloc[0],
                        gap_status["seq"].iloc[0], gap_status["align_status"].iloc[0], last_label,
                        gap_status["gap_length"].iloc[0]
                    )
                )

            else:
                summary.append(
                (
                    gap_id, gap_status["sample"].iloc[0], gap_status["asm_unit"].iloc[0],
                    gap_status["seq"].iloc[0], gap_status["align_status"].iloc[0], "not_closed",
                    gap_status["gap_length"].iloc[0]
                )
            )

    summary = pd.DataFrame.from_records(
        summary, columns=["gap_id", "sample", "asm_unit", "asm_seq", "aln_status", "max_stringency", "gap_length"]
    )

    return summary


def load_assembly_karyotype(karyotype, sample, assembly):

    norm_karyo = {
        "f": "female",
        "female": "female",
        "m": "male",
        "male": "male",
        "any": "any",
        "unknown": "any",
        "undefined": "any"
    }
    try:
        assm_karyotype = norm_karyo[karyotype.lower()]
    except KeyError:
        karyotype_file = pl.Path(karyotype).resolve(strict=True)
        karyotypes = pd.read_csv(karyotype_file, sep="\t", header=0, comment="#")
        get_sample = karyotypes["sample"] == sample
        get_assm = karyotypes["asm_unit"] == assembly
        get_karyo = get_sample & get_assm
        if not get_karyo.any():
            sys.stderr.write(f"\nWARNING:\nNo karyotype found for {sample} / {assembly}, using 'any'...\n\n")
            assm_karyotype = "any"
        else:
            assm_karyotype = karyotypes.loc[get_karyo, "karyotype"].iloc[0]
            assert assm_karyotype in norm_karyo
    return assm_karyotype


def main():

    args = parse_command_line()

    karyotype = load_assembly_karyotype(args.karyotype, args.sample, args.assembly)

    gaps = load_gap_ids(args.gap_file, karyotype)
    chrom_assign = load_chrom_assign(args.chrom_assign)

    gap_table = process_gap_alnblock_table(
        args.sample,
        args.assembly,
        args.block_intersect,
        args.flagger_regions,
        args.nucfreq_regions,
        gaps,
        chrom_assign
    )

    summary = summarize_status(gap_table)
    summary = gaps.merge(summary, on="gap_id", how="outer")
    chrom_count_stats = summary.groupby(["chrom", "aln_status", "max_stringency"]).size().reset_index()
    chrom_count_stats.rename({0: "count"}, axis=1, inplace=True)
    chrom_size_stats = summary.groupby(["chrom", "aln_status", "max_stringency"])["gap_length"].median().reset_index()
    chrom_size_stats.rename({"gap_length": "median_gap_length"}, axis=1, inplace=True)
    chrom_size_stats["median_gap_length"] = chrom_size_stats["median_gap_length"].round(0).astype(int)
    chrom_count_stats = chrom_count_stats.merge(chrom_size_stats, on=["chrom", "aln_status", "max_stringency"], how="outer")

    genome_count_stats = summary.groupby(["aln_status", "max_stringency"]).size().reset_index()
    genome_count_stats.rename({0: "count"}, axis=1, inplace=True)
    genome_size_stats = summary.groupby(["aln_status", "max_stringency"])["gap_length"].median().reset_index()
    genome_size_stats.rename({"gap_length": "median_gap_length"}, axis=1, inplace=True)
    genome_size_stats["median_gap_length"] = genome_size_stats["median_gap_length"].round(0).astype(int)
    genome_count_stats = genome_count_stats.merge(genome_size_stats, on=["aln_status", "max_stringency"], how="outer")
    genome_count_stats["chrom"] = "genome"

    genome_stats = pd.concat([genome_count_stats, chrom_count_stats], axis=0, ignore_index=False)
    genome_stats["sample"] = args.sample
    genome_stats["asm_unit"] = args.assembly
    genome_stats["sex"] = karyotype
    genome_stats.sort_values(["chrom", "aln_status", "max_stringency"], inplace=True)
    genome_stats = genome_stats[[
        "sample", "asm_unit", "sex",
        "chrom", "aln_status", "max_stringency",
        "count", "median_gap_length"
    ]].copy()

    if args.out_summary is not None:
        args.out_summary.parent.mkdir(exist_ok=True, parents=True)
        genome_stats.to_csv(args.out_summary, sep="\t", header=True, index=False)

    gaps = gaps.merge(gap_table, on="gap_id")
    gaps["sex"] = karyotype
    gaps.sort_values(["chrom", "gap_start", "gap_end"], inplace=True)
    gaps = gaps[[
        "chrom", "gap_start", "gap_end", "gap_id", "gap_length",
        "sample", "asm_unit", "sex", "seq", "start", "end", "overlap_pct",
        "align_status", "flagger_bp", "flagger_pct", "nucfreq_bp", "nucfreq_pct"
    ]].copy()

    if args.out_table is not None:
        args.out_table.parent.mkdir(exist_ok=True, parents=True)
        gaps.to_csv(args.out_table, sep="\t", header=True, index=False)

    return 0


if __name__ == "__main__":
    main()
