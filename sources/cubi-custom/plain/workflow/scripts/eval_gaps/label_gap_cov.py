#!/usr/bin/env python3

import argparse as argp
import pathlib as pl

import pandas as pd


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--input-table", "-it",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="input_table",
    )

    parser.add_argument(
        "--aln-label", "-l",
        type=str,
        choices=["bplvl", "span"],
        dest="aln_label",
        required=True
    )

    parser.add_argument(
        "--output-table", "-ot",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="output_table"
    )

    args = parser.parse_args()

    return args


def get_stage_two_labels():

    relabels = {
        ('hap:closed', 'un:ignore'): "closed|closed",
        ('hap:closed',): "closed|closed",
        ('hap:closed', 'hap:open', 'un:ignore'): "closed|open",
        ('hap:closed', 'hap:open'): "closed|open",
        ('hap:closed', 'hap:covered', 'un:ignore'): "closed|covered",
        ('hap:closed', 'hap:covered'): "closed|covered",
        ('hap:covered', 'un:ignore'): "covered|covered",
        ('hap:covered',): "covered|covered",
        ('hap:covered', 'hap:open', 'un:ignore'): "asm_err|asm_err",
        ('hap:covered', 'hap:open'): "asm_err|asm_err",
        ('hap:closed', 'hap:open', 'un:phase_hap'): "closed|ps_err",
        ('hap:closed', 'hap:open'): "closed|open",
        ('hap:open', 'un:ignore'): "open|open",
        ('hap:open',): "open|open",
        ('hap:closed', 'un:phase_hap'): "asm_err|asm_err",
        ('hap:open', 'un:error'): "ps_err|ps_err",
        ('hap:open',): "open|open",
        ('hap:closed', 'un:error'): "asm_err|asm_err",
        ('hap:closed',): "closed|closed",
        ('hap:closed', 'hap:open', 'un:phase_dip'): "ps_err|ps_err",
        ('hap:closed', 'hap:open'): "closed|open",
        ('hap:closed', 'hap:covered', 'un:error'): "asm_err|asm_err",
        ('hap:closed', 'hap:covered'): "closed|covered",
        ('hap:open', 'un:covered'): "ps_err|ps_err",
        ('hap:open',): "open|open",
        ('hap:closed', 'hap:covered', 'un:covered'): "asm_err|asm_err",
        ('hap:closed', 'hap:covered'): "closed|covered",
        ('hap:covered', 'hap:open', 'un:covered'): "asm_err|asm_err",
        ('hap:covered', 'hap:open'): "ps_err|ps_err",
        ('hap:open', 'un:phase_dip'): "ps_err|ps_err",
        ('hap:open',): "open|open",
        ('hap:closed', 'un:phase_dip'): "asm_err|asm_err",
        ('hap:closed',): "closed|closed",
        ('hap:covered', 'hap:open', 'un:error'): "asm_err|asm_err",
        ('hap:covered', 'hap:open'): "covered|open",
        ('hap:closed', 'hap:open', 'un:error'): "asm_err|asm_err",
        ('hap:closed', 'hap:open'): "closed|open",
        ('hap:closed', 'un:covered'): "asm_err|asm_err",
        ('hap:closed',): "closed|closed",
        ('hap:closed', 'hap:covered', 'un:phase_hap'): "closed|asm_err",
        ('hap:closed', 'hap:covered'): "closed|covered",
        ('hap:open', 'un:phase_hap'): "asm_err|asm_err",
        ('hap:open', ): "open|open",
        ('hap:covered', 'un:phase_hap'): "asm_err|asm_err",
        ('hap:covered',): "covered|covered",
        ('hap:covered', 'hap:open', 'un:phase_hap'): "covered|ps_err",
        ('hap:covered', 'hap:open'): "covered|open",
        ('hap:covered', 'un:covered'): "asm_err|asm_err",
        ('hap:covered',): "covered|covered",
        ('hap:closed', 'hap:open', 'un:covered'): "closed|ps_err",
        ('hap:closed', 'hap:open'): "closed|open",
        ('hap:covered', 'un:error'): "asm_err|asm_err",
        ('hap:covered',): "covered|covered",
        ('hap:covered', 'hap:open', 'un:phase_dip'): "covered|ps_err",
        ('hap:covered', 'hap:open'): "covered|open",
        ('hap:closed', 'hap:covered', 'un:phase_dip'): "asm_err|asm_err",
        ('hap:closed', 'hap:covered'): "closed|covered",
        ('hap:covered', 'un:phase_dip'): "asm_err|asm_err",
        ('hap:covered',): "covered|covered"
    }

    return relabels


def stage_one_labeling(table, aln_label):

    complex_labels = []
    for row in table.itertuples():
        au_label = "hap" if "hap" in row.asm_unit else "un"

        if row.cov_one == row.num_windows:
            label = "closed"
            if au_label == "un":
                label = "phase_hap"
        elif row.cov_two == row.num_windows and au_label == "un":
            label = "phase_dip"
        elif row.cov_tnz == row.num_windows:
            label = "covered"
        elif row.cov_tnz != row.num_windows:
            label = "open"
            if au_label == "un":
                if row.cov_tnz == 0:
                    label = "ignore"
                else:
                    label = "error"
        else:
            label = "other"

        complex_labels.append(f"{au_label}:{label}")

    table[f"gap_state_complex_{aln_label}"] = complex_labels

    return table


def stage_two_labeling(table, aln_label):

    stage2_labels = get_stage_two_labels()

    complex_label = f"gap_state_complex_{aln_label}"
    assert complex_label in table

    simple_label_column = f"gap_state_simple_{aln_label}"
    simple_labels = []

    for (sample, gap_id), cov_infos in table.groupby(["sample", "gap_id"]):
        indices = cov_infos.index
        label = tuple(sorted(set(cov_infos[complex_label])))
        simple_label = stage2_labels[label]
        simple_labels.extend(
            [(i, simple_label) for i in indices]
        )

    simple_labels = pd.DataFrame.from_records(
        simple_labels, columns=["index", simple_label_column]
    )
    simple_labels.set_index("index", inplace=True)
    table = table.merge(simple_labels, left_index=True, right_index=True)

    return table


def main():

    args = parse_command_line()

    table = pd.read_csv(args.input_table, sep="\t", header=0)
    table = stage_one_labeling(table, args.aln_label)

    table = stage_two_labeling(table, args.aln_label)

    args.output_table.parent.mkdir(exist_ok=True, parents=True)
    table.to_csv(args.output_table, sep="\t", header=True, index=False)

    return 0


if __name__ == "__main__":
    main()
