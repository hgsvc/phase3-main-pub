#!/usr/bin/env python3

import argparse as argp
import contextlib as ctl
import io
import pathlib as pl
import sys
import time

import dnaio
import xopen


def parse_command_line():

    parser = argp.ArgumentParser(prog="seq_hpc.py")

    parser.add_argument(
        "--input",
        "-in",
        "-i",
        default="stdin",
        type=str,
        dest="input",
        help="Full path to FASTA/FASTQ input file or 'stdin'. Default: stdin"
    )

    parser.add_argument(
        "--in-format",
        "-fmt",
        default="fasta",
        type=str,
        choices=["fasta", "fastq"],
        dest="input_format",
    )

    parser.add_argument(
        "--force-format",
        "-f",
        action="store_true",
        default=False,
        dest="force_format",
        help="Only for file input: force format, otherwise use auto-detection."
    )

    parser.add_argument(
        "--output",
        "-out",
        "-o",
        default="stdout",
        type=str,
        dest="output",
        help="Full path to FASTA output file or 'stdout'. Default: stdout"
    )

    parser.add_argument(
        "--cmap-table",
        "-ct",
        default=None,
        type=lambda x: pl.Path(x).resolve(),
        dest="cmap_table",
    )

    parser.add_argument(
        "--report",
        "-r",
        action="store_true",
        default=False,
        dest="report",
        help="If set, write brief report to sys.stderr. Default: False"
    )

    parser.add_argument(
        "--buffer-size",
        "-b",
        type=int,
        default=0,
        dest="buffer_size",
        help="Set buffer size (#char) or 0 for no buffer. Default: 0"
    )

    parser.add_argument(
        "--skip-self-test",
        action="store_true",
        default=False,
        dest="skip_self_test",
        help="Skip self test of implementation at startup. Default: False"
    )

    parser.add_argument(
        "--break-after",
        "-ba",
        type=int,
        default=0,
        dest="break_after",
        help="Abort operations after this many sequences; 0 to process all. Default: 0"
    )

    args = parser.parse_args()

    return args


def homopolymer_compress(sequence):
    """NB: this would fail on empty input

    Args:
        sequence (str): sequence to hpc
    Returns:
        str: hpc sequence
    """
    uncompressed_length = len(sequence)
    hpc_seq = "".join(a if a != b else "" for a,b in zip(sequence[:-1], sequence[1:]))
    if not hpc_seq:
        hpc_seq = sequence[0]
    if hpc_seq[-1] != sequence[-1]:
        hpc_seq += sequence[-1]
    compressed_length = len(hpc_seq)
    ratio = round(compressed_length / uncompressed_length, 3)
    assert ratio <= 1
    # None in return for uniform interface
    # with cmap version
    return hpc_seq, ratio, None


def homopolymer_compress_with_cmap(sequence):
    """NB: this would fail on empty input

    This version is about ~50% slower, but creates
    a 1-to-1 coordinate map between sequences
    on the fly.

    Args:
        sequence (str): sequence to hpc
    Returns:
        str: hpc sequence
    """
    uncompressed_length = len(sequence)
    hpc_seq = sequence[0]
    cmap = []
    raw_block_start = 0
    raw_block_end = 1
    hpc_block_start = 0
    hpc_block_end = 1
    hpc_stretch = False
    for char in sequence[1:]:
        if char == hpc_seq[-1]:
            raw_block_end += 1
            hpc_stretch = True
        elif not hpc_stretch:
            hpc_seq += char
            raw_block_end += 1
            hpc_block_end += 1
        else:
            # homopolymer stretch ended
            cmap.append((hpc_block_start, hpc_block_end, raw_block_start, raw_block_end))
            hpc_seq += char
            raw_block_start = raw_block_end
            raw_block_end += 1
            hpc_block_start = hpc_block_end
            hpc_block_end += 1
            hpc_stretch = False
    cmap.append((hpc_block_start, hpc_block_end, raw_block_start, raw_block_end))
    compressed_length = len(hpc_seq)
    ratio = round(compressed_length / uncompressed_length, 3)

    if not cmap:
        raise ValueError("Failed to generate coordinate map")

    assert ratio <= 1
    return hpc_seq, ratio, cmap


def check_implementation():

    test_strings = [
        "AAA", "ACGT", "ACAA", "XXXNNX",
        "AACGT", "ACCCCGGGGT","ACGTTT"
    ]

    result_strings = [
        "A", "ACGT", "ACA", "XNX",
        "ACGT", "ACGT", "ACGT"
    ]

    assert all(
        [homopolymer_compress(test)[0] == result]
        for test, result in zip(test_strings, result_strings)
    )

    assert all(
        [homopolymer_compress_with_cmap(test)[0] == result]
        for test, result in zip(test_strings, result_strings)
    )

    return


def main():

    args = parse_command_line()

    if not args.skip_self_test:
        _ = check_implementation()

    if args.input not in ["stdin", "-", "/dev/stdin", ""]:
        infile = pl.Path(args.input).resolve(strict=True)
    else:
        infile = sys.stdin.buffer

    if not args.force_format:
        input_format = None
    else:
        input_format = args.input_format

    if args.output not in ["stdout", "-", "/dev/stdout", ""]:
        outfile = pl.Path(args.output).resolve()
        outfile.parent.mkdir(exist_ok=True, parents=True)
    else:
        outfile = sys.stdout.buffer


    if args.cmap_table is not None:
        hpc = homopolymer_compress_with_cmap
        args.cmap_table.parent.mkdir(exist_ok=True, parents=True)
        coordinate_map = io.StringIO()
        cmap_table_header = "\t".join(
            ["#name", "hpc_start", "hpc_end", "seqnum", "plain_start", "plain_end", "hpc_ratio"]
        ) + "\n"
        coordinate_map.write(cmap_table_header)
    else:
        hpc = homopolymer_compress
        coordinate_map = None

    buffer_limit = args.buffer_size
    buffered = 0
    use_buffer = buffer_limit > 0
    buffer_obj = None
    write_buffer = None

    avg_hpc_ratio = 0.
    process_start = time.perf_counter()

    with ctl.ExitStack() as exs:
        read_input = exs.enter_context(
            dnaio.open(infile, fileformat=input_format, mode="r")
        )

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

        if args.cmap_table is not None:
            cmap_table = exs.enter_context(xopen.xopen(args.cmap_table, mode="w"))

        for record_num, seq_record in enumerate(read_input, start=1):
            hpc_seq, ratio, cmap = hpc(seq_record.sequence.upper())
            avg_hpc_ratio += ratio
            if cmap is not None:
                [
                    coordinate_map.write(
                        (
                            f"{seq_record.name}\t"
                            f"{cmap_block[0]}\t"
                            f"{cmap_block[1]}\t"
                            f"{record_num}\t"
                            f"{cmap_block[2]}\t"
                            f"{cmap_block[3]}\t"
                            f"{ratio}\n"
                        )
                    )
                    for cmap_block in cmap
                ]


            if use_buffer:
                write_buffer.write(seq_record.name, hpc_seq)
                buffered += len(seq_record.name)
                buffered += len(hpc_seq)
                if buffered > buffer_limit:
                    write_output.write(buffer_obj.getvalue())
                    buffer_obj = io.BytesIO()
                    write_buffer = dnaio.open(buffer_obj, mode="w", fileformat="fasta")
                    buffered = 0
                    if cmap is not None:
                        cmap_table.write(coordinate_map.getvalue())
                        coordinate_map = io.StringIO()

            else:
                write_output.write(seq_record.name, hpc_seq)
                if cmap is not None:
                    cmap_table.write(coordinate_map.getvalue())
                    coordinate_map = io.StringIO()

            if args.break_after > 0 and record_num > args.break_after:
                break


        if buffered > 0:
            write_output.write(buffer_obj.getvalue())
            if cmap is not None:
                cmap_table.write(coordinate_map.getvalue())

    process_end = time.perf_counter()
    total_time = round(process_end - process_start, 3)

    avg_hpc_ratio = round(avg_hpc_ratio / record_num, 3)
    if args.report:
        sys.stderr.write(
            f"\n=== seq_hpc.py report ==="
            f"\nProcessed records: {record_num}"
            f"\nAverage homopolymer-compression ratio: {avg_hpc_ratio}"
            f"\nTotal processing time: ~{total_time} sec\n\n"
        )

    return 0


if __name__ == "__main__":
    main()
