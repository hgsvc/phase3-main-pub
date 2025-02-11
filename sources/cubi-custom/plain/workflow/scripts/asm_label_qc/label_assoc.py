#!/usr/bin/env python3

import argparse as argp
import decimal
import pathlib as pl

import pandas as pd
import scipy.stats as stats
import scipy.stats.contingency as contin
import numpy as np


def parse_cli():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--augmented-regions", "-ar",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="regions",
        help="Error flagged regions augmented with annotations"
    )

    parser.add_argument(
        "--output", "-o",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="output",
        help="Output table (TSV)"
    )

    args = parser.parse_args()
    return args


def extract_flags(labels):
    flags = set(l.split("::")[0] for l in labels.split(","))
    return flags


def main():

    args = parse_cli()
    table = pd.read_csv(args.regions, sep="\t", header=0)

    # drop seq/contigs w/o any annotated errors (mostly very short / garbage)
    table = table.loc[~table["labels"].str.contains("no-labels"), :].copy()

    error_flags = set().union(*table["labels"].apply(extract_flags))
    annotations = ["telomere", "centromere", "segdup98"]
    assert all(f"overlaps_{a}" in table.columns for a in annotations)

    stat_rows = []
    raw_pvalues = []
    for error_flag in error_flags:
        for ann in annotations:
            has_flag = table["labels"].str.contains(error_flag)
            overlaps = table[f"overlaps_{ann}"] == 1

            # add fudge factor of 1 to avoid zero div issues
            flag_overlaps = (has_flag & overlaps).sum() + 1
            flag_not_overlaps = (has_flag & ~overlaps).sum() + 1
            other_overlaps = (~has_flag & overlaps).sum() + 1
            other_not_overlaps = (~has_flag & ~overlaps).sum() + 1

            # 2x2 table
            # cases [overlaps annotation]: w/ err label /// w/o err label
            # non-cases [no overlap]: w/ label /// w/o label
            contingency_table = np.array([[flag_overlaps, other_overlaps], [flag_not_overlaps, other_not_overlaps]], dtype=int)

            fet_or, fet_pv = stats.fisher_exact(contingency_table, alternative="two-sided")
            or_result = contin.odds_ratio(contingency_table, kind="sample")
            odds_ratio = or_result.statistic
            assert np.isclose(fet_or, odds_ratio, atol=1e-5)
            ci_low, ci_high = or_result.confidence_interval(0.95, alternative="two-sided")
            if np.isnan(ci_low) or np.isnan(ci_high):
                # should not happen b/c of fudge count of 1
                raise ValueError(error_flag, ann)

            stat_rows.append(
                [error_flag, ann, round(odds_ratio, 3), round(ci_low, 3), round(ci_high, 3), f"{decimal.Decimal(fet_pv):.2E}"]
            )
            raw_pvalues.append(fet_pv)

    raw_pvalues = np.array(raw_pvalues, dtype=float)
    adj_pvalues = stats.false_discovery_control(raw_pvalues, method="by")
    alpha = 1e-3
    for pos, pv in enumerate(adj_pvalues):
        if pv < alpha:
            label = "sig."
        else:
            label = "n.s."
        fmt_pv = f"{decimal.Decimal(pv):.2E}"
        stat_rows[pos].extend([fmt_pv, label, alpha])

    df = pd.DataFrame.from_records(
        stat_rows,
        columns=["error_label", "annotation", "odds_ratio", "or_ci_low", "or_ci_high", "fet_pvalue", "adj_pvalue", "sig_label", "alpha"]
    )
    df.sort_values(["error_label", "annotation"], inplace=True)
    args.output.parent.mkdir(exist_ok=True, parents=True)
    df.to_csv(args.output, sep="\t", header=True, index=False)

    return 0


if __name__ == "__main__":
    main()
