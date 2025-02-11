#!/usr/bin/env python3

import argparse as argp
import json
import pathlib as pl
import shutil


KNOWN_OUTPUT_FILES = [
    (
        "haplotype1.fasta",
        "hap1_fasta",
        "Haplotype1 assembly."
    ),
    (
        "haplotype2.fasta",
        "hap2_fasta",
        "Haplotype2 assembly."
    ),
    (
        "unassigned.fasta",
        "unassigned_fasta",
        "Unassigned assembly contigs (unphased contigs)."
    ),
    (
        "disconnected.fasta",
        "disconnected_fasta",
        "Short contigs not connected to larger components, likely garbage.",
    ),
    ("ebv.exemplar.fasta", "ebv_repr", "Representative EBV sequence."),
    ("ebv.fasta", "ebv_fasta", "All sequences identified as EBV."),
    ("fasta", "wg_fasta", "Whole-genome assembly."),
    ("hifi-coverage.csv", "hifi_cov", "HiFi coverage estimates on unitig level."),
    (
        "homopolymer-compressed.gfa",
        "wg_gfa_hpc",
        "Whole-genome assembly graph with sequences in HPC space",
    ),
    (
        "homopolymer-compressed.layout",
        "wg_layout",
        "Whole-genome graph layout information.",
    ),
    (
        "homopolymer-compressed.noseq.gfa",
        "wg_gfa_noseq",
        "Whole-genome assembly graph with sequences in HPC space",
    ),
    ("mito.exemplar.fasta", "mito_repr", "Representative mitochondrial sequence."),
    ("mito.fasta", "mito_fasta", "All sequences identified as mitochondrial."),
    ("ont-coverage.csv", "ont_cov", "ONT coverage estimates on unitig level."),
    ("rdna.exemplar.fasta", "rdna_repr", "Representative rDNA sequence."),
    ("rdna.fasta", "rdna_fasta", "All sequences identified as rDNA."),
    ("scfmap", "scf_map", "Scaffold mapping."),
    ("colors.csv", "node_coloring", "Rukki node coloring (only in trio mode)."),
    ("paths.tsv", "rukki_paths", "Rukki paths file (only in trio mode)."),

]

KNOWN_OUTPUT_FILES = dict((t[0], t[1:]) for t in KNOWN_OUTPUT_FILES)


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--verkko-wd",
        "-w",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="work_dir",
        required=True,
        help="Path to Verkko working directory (top-level).",
    )
    parser.add_argument(
        "--base-dir",
        "-b",
        required=True,
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="base_dir",
    )
    parser.add_argument(
        "--sample", "-s", type=str, default=None, help="Add sample name to output file"
    )
    parser.add_argument(
        "--delete-logs",
        "-d",
        action="store_true",
        default=False,
        dest="delete_logs",
        help="Delete the subfolder log/ or logs/ if it exists in the working directory.",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=lambda x: pl.Path(x),
        dest="output",
        help="Path to output json file storing file paths.",
    )
    args = parser.parse_args()

    return args


def main():

    args = parse_command_line()

    file_collection = dict()
    if args.sample is not None:
        file_collection["sample"] = args.sample
    for output_file in args.work_dir.glob("assembly.*"):
        file_ext = output_file.name.split(".", 1)[-1]
        file_key, file_desc = KNOWN_OUTPUT_FILES[file_ext]
        file_abs_path = output_file.resolve(strict=True)
        file_rel_path = file_abs_path.relative_to(args.base_dir)
        file_info = {
            "desc": file_desc,
            "abs_path": str(file_abs_path),
            "rel_path": str(file_rel_path),
        }
        file_collection[file_key] = file_info

    if "wg_fasta" not in file_collection:
        raise RuntimeError(f"Incomplete Verkko run detected: {args.work_dir}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as dump:
        json.dump(file_collection, dump, ensure_ascii=True, indent=2)

    if args.delete_logs:
        check_log_folders = ["log", "logs"]
        for check_log in check_log_folders:
            potential_log_path = args.work_dir.joinpath(check_log)
            if potential_log_path.is_dir():
                shutil.rmtree(potential_log_path, ignore_errors=True)

    return 0


if __name__ == "__main__":
    main()
