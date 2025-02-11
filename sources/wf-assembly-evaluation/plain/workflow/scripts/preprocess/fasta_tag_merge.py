#!/usr/bin/env python3

import argparse as argp
import collections as col
import contextlib as ctl
import io
import pathlib as pl
import re
import sys
import time

import xopen
import dnaio


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--input-files",
        "--input",
        "-i",
        nargs="+",
        type=lambda x: pl.Path(x).resolve(strict=True),
        required=True,
        dest="input_files",
        help="Full path(s) to FASTA input files."
    )

    parser.add_argument(
        "--seq-tags",
        "--tags",
        "-t",
        type=lambda x: pl.Path(x).resolve(strict=True),
        required=False,
        dest="seq_tags",
        default=None,
        help="Path to sequence tags file (TSV: filename <TAB> tag). Default: <NONE>"
    )

    parser.add_argument(
        "--exclude", "--ignore",
        "--skip", "--discard",
        type=lambda x: pl.Path(x).resolve(),
        default=None,
        dest="skip_seqs",
        help="Text file listing sequences (by name) to be ignored / skipped over. Default: <NONE>"
    )

    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default="stdout"
    )

    parser.add_argument(
        "--buffer-size",
        "-b",
        type=int,
        default=0,
        dest="buffer_size"
    )

    # 2024-01-16
    # the following argument was introduced
    # to make the NCBI FCS adaptor screening work
    # again for some assemblies that contained
    # sequences shorter than 10 bp (origin unknown)
    parser.add_argument(
        "--skip-scraps",
        "-s",
        type=int,
        default=10,
        dest="skip_scraps",
        help=(
            "Skip over scraps (very short sequence fragments) "
            "shorter than N bp in the input. Default: 10 (bp)"
        )
    )

    parser.add_argument(
        "--report",
        "-r",
        action="store_true",
        default=False,
        dest="report",
        help="If set, write brief report to sys.stderr. Default: False"
    )

    args = parser.parse_args()

    return args


def read_sequence_tags(file_path):

    if file_path is None:
        tags = col.defaultdict(str)
    else:
        valid_tag = re.compile("[A-Za-z0-9]+")
        tags = dict()
        seen_tags = set()
        with open(file_path, "r") as table:
            for line in table:
                filename, tag = line.strip().split()
                assert filename not in tags
                if valid_tag.match(tag) is None:
                    raise ValueError(
                        f"Invalid sequence tag: {tag}\n"
                        "(Allowed: A-Z, a-z, 0-9)"
                    )
                assert tag not in seen_tags
                tags[filename] = f".{tag}"
                seen_tags.add(tag)
    return tags


def load_discard_seq_names(file_path):

    if file_path is None:
        skip_seqs = set()
    else:
        file_path = file_path.resolve(strict=True)
        with xopen.xopen(file_path) as listing:
            skip_seqs = set(listing.read().strip().split())
        if len(skip_seqs) < 1:
            raise ValueError(
                f"No sequence names to skip loaded from file: {file_path}"
            )
    return skip_seqs


def main():

    args = parse_command_line()

    use_buffer = args.buffer_size > 0
    buffer_limit = args.buffer_size
    buffer_obj = None
    write_buffer = None
    buffered = 0

    file_tags = read_sequence_tags(args.seq_tags)

    skip_seqs = load_discard_seq_names(args.skip_seqs)

    if args.output in ["stdout", "-", "/dev/stdout", ""]:
        outfile = sys.stdout.buffer
    else:
        outfile = pl.Path(args.output).resolve()
        outfile.parent.mkdir(exist_ok=True, parents=True)

    processed_records = 0
    skipped_scraps = 0  # see above in arg parser
    skipped_seqs = 0
    process_start = time.perf_counter()
    check_uniq_seqnames = col.defaultdict(list)
    with ctl.ExitStack() as exs:

        if use_buffer:
            buffer_obj = io.BytesIO()
            write_buffer = dnaio.open(buffer_obj, mode="w", fileformat="fasta")
            write_output = exs.enter_context(
                xopen.xopen(outfile, mode="wb", compresslevel=5)
            )
        else:
            write_output = exs.enter_context(
                dnaio.open(
                    outfile, fileformat="fasta", mode="w",
                    qualities=False, compression_level=5
                )
            )

        for input_file in args.input_files:
            file_tag = file_tags[input_file.name]
            with dnaio.open(input_file, mode="r") as fasta:
                for record in fasta:
                    processed_records += 1
                    if len(record.sequence) < args.skip_scraps:
                        skipped_scraps += 1
                        continue
                    if record.name in skip_seqs:
                        skipped_seqs += 1
                        continue
                    check_uniq_seqnames[record.name].append(file_tag)
                    tagged_name = f"{record.name}{file_tag}"
                    if use_buffer:
                        write_buffer.write(tagged_name, record.sequence)
                        buffered += len(tagged_name)
                        buffered += len(record.sequence)
                        if buffered > buffer_limit:
                            write_output.write(buffer_obj.getvalue())
                            buffer_obj = io.BytesIO()
                            write_buffer = dnaio.open(buffer_obj, mode="w", fileformat="fasta")
                            buffered = 0
                    else:
                        write_output.write(tagged_name, record.sequence)

                if buffered > 0:
                    write_output.write(buffer_obj.getvalue())
                    buffer_obj = io.BytesIO()
                    write_buffer = dnaio.open(buffer_obj, mode="w", fileformat="fasta")
                    buffered = 0

    process_middle = time.perf_counter()
    spurious_seqs = col.Counter()
    for seq, tags in check_uniq_seqnames.items():
        if len(tags) > 1:
            spurious_seqs[tuple(sorted(tags))] += 1

    process_end = time.perf_counter()
    tagging_time = round(process_middle - process_start, 3)
    tag_eval_time = round(process_end - process_middle, 3)
    total_time = round(process_end - process_start, 3)

    non_uniq_seq_names = "\n".join(["\t".join(["|".join(n), str(c)]) for n, c in spurious_seqs.most_common()])

    if not non_uniq_seq_names.strip():
        non_uniq_seq_names = "<NONE>"

    if args.report:
        sys.stderr.write(
            "\n\n=== fasta_tag_merge report ==="
            f"\nProcessed records: {processed_records}"
            f"\nSkipped scrap records: {skipped_scraps}"
            f"\nSkipped seq. records by name: {skipped_seqs}"
            f"\nTotal processing time: ~{total_time} sec"
            f"\nTagging time (incl. I/O): ~{tagging_time} sec"
            f"\nTag evaluation time: ~{tag_eval_time} sec"
            f"\nNon-unique sequence names were tagged:\n\n"
            f"{non_uniq_seq_names}\n\n"
        )

    return 0



if __name__ == "__main__":
    main()
