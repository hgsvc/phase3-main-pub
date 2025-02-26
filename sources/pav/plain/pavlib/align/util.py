# Alignment base operations

import numpy as np
import pandas as pd
import pysam

import pavlib.seq

# CIGAR operations
_INT_STR_SET = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'}
_CIGAR_OP_SET = {'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'}

CIGAR_M = 0
CIGAR_I = 1
CIGAR_D = 2
CIGAR_N = 3
CIGAR_S = 4
CIGAR_H = 5
CIGAR_P = 6
CIGAR_EQ = 7
CIGAR_X = 8

CIGAR_CODE_TO_CHAR = {
    CIGAR_M: 'M',
    CIGAR_I: 'I',
    CIGAR_D: 'D',
    CIGAR_N: 'N',
    CIGAR_S: 'S',
    CIGAR_H: 'H',
    CIGAR_P: 'P',
    CIGAR_EQ: '=',
    CIGAR_X: 'X'
}

# Indices for tuples returned by trace_cigar_to_zero()
TC_INDEX = 0
TC_OP_LEN = 1
TC_OP_CODE = 2
TC_DIFF_CUM = 3
TC_DIFF = 4
TC_EVENT_CUM = 5
TC_EVENT = 6
TC_SUB_BP = 7
TC_QRY_BP = 8
TC_CLIPS_BP = 9
TC_CLIPH_BP = 10

TRIM_DESC = {
    'none': 'No trimming',
    'qry': 'Query-only trimming',
    'qryref': 'Query-Reference trimming'
}


def align_bed_to_depth_bed(df, df_fai=None):
    """
    Get a BED file of alignment depth from an alignment BED.

    Output columns:
    * #CHROM: Reference chromosome.
    * POS: Reference start position (BED coordinates).
    * END: Reference end position (BED coordinates).
    * DEPTH: Number of alignment records. Integer, may be 0.
    * QUERY: Query names. Comma-separated list if multiple queries.
    * INDEX: Query indexes in the same order as QUERIES (corresponding INDEX column in df).

    If `df_fai` is not 'None', then depth records extend from 0 to the end of the chromosome. If the first record does
    not reach position 0, a 0-depth record is added over that region. Similarly, if the last record does not reach the
    end of the chromosome, a 0-depth record is added over that region. If `df_fai` is `None`, then no padding is done
    to the start or end of the chromosome and 0-depth records only appear between alignment records.

    :param df: Alignment BED file.
    :param df_fai: FAI series (keys = reference sequence names, values = sequence lengths).

    :return: A Pandas DataTable with depth across all reference loci.
    """

    # Get a list of alignment events (start and end alignment record)

    # Build an ordered list of alignments.
    # align_list is a list of tuples:
    #  0) Chromosome
    #  1) Event position
    #  2) Event type (1 is align start, 0 is align end)
    #  3) Align record index. Removes the correct alignment from the aligned query list when an alignment record ends.
    #  4) Query ID: List of query IDs in the aligned region (comma-separated list). Each alignment start goes to the
    #     end of the list, and each alignment end removes the element from the list that it added even if the query
    #     ID is the same.

    align_list = list()

    for index, row in df.iterrows():
        align_list.append(
            (str(row['#CHROM']), row['POS'], 1, row['INDEX'], row['QRY_ID'], row['INDEX'])
        )

        align_list.append(
            (str(row['#CHROM']), row['END'], 0, row['INDEX'], row['QRY_ID'], row['INDEX'])
        )

    if not align_list:
        raise RuntimeError('No alignments to process')

    align_list = sorted(align_list)

    if align_list[0][2] != 1:
        raise RuntimeError(f'First alignment is not a record start: {",".join([str(val) for val in align_list[0]])}')

    # Setup to process BED records
    df_bed_list = list()

    last_chrom = None
    last_pos = None
    qry_list = list()

    # Write empty records for reference sequences with no alignments
    if df_fai is not None:
        seq_fai_list = list(sorted(df_fai.index.astype(str)))

        while len(seq_fai_list) > 0 and seq_fai_list[0] < align_list[0][0]:
            df_bed_list.append(
                pd.Series(
                    [seq_fai_list[0], 0, df_fai[seq_fai_list[0]], 0, '', ''],
                    index=['#CHROM', 'POS', 'END', 'DEPTH', 'QRY_ID', 'INDEX']
                )
            )

            seq_fai_list = seq_fai_list[1:]

        if len(seq_fai_list) == 0 or seq_fai_list[0] != align_list[0][0]:
            raise RuntimeError(f'Missing {align_list[0][0]} in FAI or out of order')

    else:
        seq_fai_list = None

    # Process BED records
    for chrom, pos, event, index, qry, row_index in align_list:

        # Chromosome change
        if chrom != last_chrom:

            # Check sanity
            if qry_list:
                raise RuntimeError(f'Switched chromosome ({last_chrom} > {chrom}) with open queries: {", ".join(qry_list)}')

            # Check chromosome order
            if df_fai is not None:

                # Check chromosome in FAI
                if chrom not in df_fai.index:
                    raise RuntimeError(f'Missing chromosome in reference FAI index: {chrom}')

                # Write empty records for reference sequences with no alignments
                while len(seq_fai_list) > 0 and seq_fai_list[0] < chrom:
                    df_bed_list.append(
                        pd.Series(
                            [seq_fai_list[0], 0, df_fai[seq_fai_list[0]], 0, '', ''],
                            index=['#CHROM', 'POS', 'END', 'DEPTH', 'QRY_ID', 'INDEX']
                        )
                    )

                    seq_fai_list = seq_fai_list[1:]

                if len(seq_fai_list) == 0 or seq_fai_list[0] != chrom:
                    raise RuntimeError(f'Missing {chrom} in FAI or out of order')

                seq_fai_list = seq_fai_list[1:]

                # Add record up to end of chromosome
                if last_chrom is not None:
                    if last_pos > df_fai[last_chrom]:
                        raise RuntimeError(f'Last END position for chromosome {last_chrom} is greater than chromosome length: {last_pos} > {df_fai[last_chrom]}')

                    if last_pos < df_fai[last_chrom]:
                        df_bed_list.append(
                            pd.Series(
                                [
                                    last_chrom, last_pos, df_fai[last_chrom], len(qry_list),
                                    ','.join([val[1] for val in qry_list]),
                                    ','.join([str(val[0]) for val in qry_list])
                                ],
                                index=['#CHROM', 'POS', 'END', 'DEPTH', 'QRY_ID', 'INDEX']
                            )
                        )

            # Add chrom:0-pos record
            if df_fai is not None and pos > 0:
                df_bed_list.append(
                    pd.Series(
                        [
                            chrom, 0, pos, len(qry_list),
                            ','.join([val[1] for val in qry_list]),
                            ','.join([str(val[0]) for val in qry_list])
                        ],
                        index=['#CHROM', 'POS', 'END', 'DEPTH', 'QRY_ID', 'INDEX']
                    )
                )

            # Set state for the next alignment
            last_chrom = chrom
            last_pos = pos

        # Check position sanity
        if last_pos > pos:
            raise RuntimeError(f'Position out of order: {last_chrom}:{last_pos} > {chrom}:{pos}')

        # Write last record
        if pos > last_pos:
            df_bed_list.append(
                pd.Series(
                    [
                        chrom, last_pos, pos, len(qry_list),
                        ','.join([val[1] for val in qry_list]),
                        ','.join([str(val[0]) for val in qry_list])
                    ],
                    index=['#CHROM', 'POS', 'END', 'DEPTH', 'QRY_ID', 'INDEX']
                )
            )

            last_pos = pos

        # Process event
        if event == 1:
            qry_list.append((index, qry))

        elif event == 0:
            n = len(qry_list)

            if n == 0:
                raise RuntimeError(f'Got END event with no queries in the list: {chrom},{pos},{event},{index},{qry}')

            qry_list = [val for val in qry_list if val != (index, qry)]

            if len(qry_list) == n:
                raise RuntimeError(f'Could not find query to END in query list: {chrom},{pos},{event},{index},{qry}')

            if len(qry_list) < n - 1:
                raise RuntimeError(f'END removed multiple queries: {chrom},{pos},{event},{index},{qry}')

        else:
            raise RuntimeError(f'Unknown event type {event}: {chrom},{pos},{event},{index},{qry}')

    # Check final state
    if qry_list:
        raise RuntimeError(f'Ended alignment records with open queries: {", ".join(qry_list)}')

    if df_fai is not None:
        if last_chrom not in df_fai.index:
            raise RuntimeError(f'Missing chromosome in reference FAI index: {last_chrom}')

        # Add final record up to end of chromosome
        if last_pos > df_fai[last_chrom]:
            raise RuntimeError(f'Last END position for chromosome {last_chrom} is greater than chromosome length: {last_pos} > {df_fai[last_chrom]}')

        if last_pos < df_fai[last_chrom]:
            df_bed_list.append(
                pd.Series(
                    [
                        last_chrom, last_pos, df_fai[last_chrom], len(qry_list),
                        ','.join([val[1] for val in qry_list]),
                        ','.join([str(val[0]) for val in qry_list])
                    ],
                    index=['#CHROM', 'POS', 'END', 'DEPTH', 'QRY_ID', 'INDEX']
                )
            )

        # Write empty records for reference sequences with no alignments
        while len(seq_fai_list) > 0:

            df_bed_list.append(
                pd.Series(
                    [seq_fai_list[0], 0, df_fai[seq_fai_list[0]], 0, '', ''],
                    index=['#CHROM', 'POS', 'END', 'DEPTH', 'QRY_ID', 'INDEX']
                )
            )

            seq_fai_list = seq_fai_list[1:]

    # Create BED file
    df_bed = pd.concat(df_bed_list, axis=1).T

    df_bed.sort_values(['#CHROM', 'POS'], inplace=True)

    return df_bed


def cigar_str_to_tuples(record):
    """
    Get an iterator for cigar operation tuples. Each tuple is (cigar-len, cigar-op).

    :param record: Alignment record or a CIGAR string.

    :return: Iterator of CIGAR operation tuples.
    """

    if type(record) == pd.Series:
        cigar = record['CIGAR']
    else:
        cigar = record

    pos = 0
    max_pos = len(cigar)

    while pos < max_pos:

        len_pos = pos

        while cigar[len_pos] in _INT_STR_SET:
            len_pos += 1

        if len_pos == pos:
            raise RuntimeError('Missing length in CIGAR string for query {} alignment starting at {}:{}: CIGAR index {}'.format(
                record['QRY_ID'], record['#CHROM'], record['POS'], pos
            ))

        if cigar[len_pos] not in _CIGAR_OP_SET:
            raise RuntimeError('Unknown CIGAR operation for query {} alignment starting at {}:{}: CIGAR operation {}'.format(
                record['QRY_ID'], record['#CHROM'], record['POS'], cigar[pos]
            ))

        yield int(cigar[pos:len_pos]), cigar[len_pos]

        pos = len_pos + 1


def match_bp(record, right_end):
    """
    Get the number of matching bases at the end of an alignment. Used by variant callers to left-align SVs through
    alignment-truncating events.

    :param record: Alignment record (from alignment BED) with CIGAR string.
    :param right_end: `True` if matching alignments from the right end of `record`, or `False` to match from
        the left end.

    :return: Minimum of the number of matched bases at the end of two alignment records.
    """

    cigar = list(cigar_str_to_tuples(record))

    if right_end:
        cigar = cigar[::-1]

    # Get match base count (CIGAR op "=") on a
    match_count = 0

    for cigar_len, cigar_op in cigar:
        if cigar_op in {4, 5}:  # Skip clipped bases: S, H
            continue

        if cigar_op == 7:  # Matched bases: =
            match_count += cigar_len

        elif cigar_op == 0:
            raise RuntimeError(
                'Detected "M" opcodes in CIGAR string for record INDEX={}: Sequence match/mismatch opcodes are required ("=", "X")'.format(
                    record['INDEX'] if 'INDEX' in record.index else '<UNKNOWN>'
                )
            )
        else:
            break  # No more bases to traverse

    return match_count


def check_record(row, df_qry_fai):
    """
    Check alignment DatFrame record for sanity. Throws exceptions if there are problems. Returns nothing if everything
    passes.

    :param row: Alignment table record (Pandas Series).
    :param df_qry_fai: Panadas Series with query names as keys and query lengths as values.
    """

    try:
        ref_bp, qry_bp, clip_h_l, clip_s_l, clip_h_r, clip_s_r = count_cigar(row)

    except Exception as ex:
        raise RuntimeError(
            'CIGAR parsing error: {} (INDEX={}, QRY={}:{}-{}, REF={}:{}-{})'.format(
                ex,
                row['INDEX'], row['QRY_ID'], row['QRY_POS'], row['QRY_END'], row['#CHROM'], row['POS'], row['END']
            )
        )

    qry_len = df_qry_fai[row['QRY_ID']]

    if 'QRY_LEN' in row.index:
        raise RuntimeError('QRY_LEN defined (was removed from PAV alignment files) (INDEX={}, QRY={}:{}-{}, REF={}:{}-{})'.format(
            row['INDEX'], row['QRY_ID'], row['QRY_POS'], row['QRY_END'], row['#CHROM'], row['POS'], row['END']
        ))

    # if row['QRY_LEN'] != qry_len:
    #     raise RuntimeError('QRY_LEN != length from FAI ({} != {}) (INDEX={}, QRY={}:{}-{}, REF={}:{}-{})'.format(
    #         row['QRY_LEN'], qry_len,
    #         row['INDEX'], row['QRY_ID'], row['QRY_POS'], row['QRY_END'], row['#CHROM'], row['POS'], row['END']
    #     ))

    # qry_map_rgn = pavlib.seq.region_from_string(row['QRY_MAP_RGN'])

    # Query and reference positions are in the right order
    if row['QRY_POS'] >= row['QRY_END']:
        raise RuntimeError('QRY_POS >= QRY_END ({} >= {}) (INDEX={}, QRY={}:{}-{}, REF={}:{}-{})'.format(
            row['QRY_POS'], row['QRY_END'],
            row['INDEX'], row['QRY_ID'], row['QRY_POS'], row['QRY_END'], row['#CHROM'], row['POS'], row['END']
        ))

    if row['POS'] >= row['END']:
        raise RuntimeError('POS >= END ({} >= {}) (INDEX={}, QRY={}:{}-{}, REF={}:{}-{})'.format(
            row['POS'], row['END'],
            row['INDEX'], row['QRY_ID'], row['QRY_POS'], row['QRY_END'], row['#CHROM'], row['POS'], row['END']
        ))

    # No negative positions
    if row['POS'] < 0:
        raise RuntimeError('POS ({}) < 0 (INDEX={}, QRY={}:{}-{}, REF={}:{}-{})'.format(
            row['POS'],
            row['INDEX'], row['QRY_ID'], row['QRY_POS'], row['QRY_END'], row['#CHROM'], row['POS'], row['END']
        ))

    if row['QRY_POS'] < 0:
        raise RuntimeError('QRY_POS ({}) < 0 (INDEX={}, QRY={}:{}-{}, REF={}:{}-{})'.format(
            row['QRY_POS'],
            row['INDEX'], row['QRY_ID'], row['QRY_POS'], row['QRY_END'], row['#CHROM'], row['POS'], row['END']
        ))

    # if qry_map_rgn.pos < 0:
    #     raise RuntimeError('QRY_MAP_RGN.pos ({}) < 0 (INDEX={}, QRY={}:{}-{}, REF={}:{}-{})'.format(
    #         qry_map_rgn.pos,
    #         row['INDEX'], row['QRY_ID'], row['QRY_POS'], row['QRY_END'], row['#CHROM'], row['POS'], row['END']
    #     ))

    # POS and END agree with length
    if row['POS'] + ref_bp != row['END']:

        raise RuntimeError(
            'END mismatch: POS + ref_bp != END ({} != {}) (INDEX={}, QRY={}:{}-{}, REF={}:{}-{})'.format(
                row['POS'] + ref_bp, row['END'],
                row['INDEX'], row['QRY_ID'], row['QRY_POS'], row['QRY_END'], row['#CHROM'], row['POS'], row['END']
            )
        )

    # Query POS and END agree with length
    if row['QRY_POS'] + qry_bp != row['QRY_END']:
        raise RuntimeError(
            'QRY_POS + qry_bp != QRY_END: {} != {} (INDEX={}, QRY={}:{}-{}, REF={}:{}-{})'.format(
                row['QRY_POS'] + qry_bp, row['QRY_END'],
                row['INDEX'], row['QRY_ID'], row['QRY_POS'], row['QRY_END'], row['#CHROM'], row['POS'], row['END']
            )
        )

    # Query ends are not longer than query lengths
    if row['QRY_END'] > qry_len:
        raise RuntimeError('QRY_END > qry_len ({} > {}) (INDEX={}, QRY={}:{}-{}, REF={}:{}-{})'.format(
            row['QRY_END'], qry_len,
            row['INDEX'], row['QRY_ID'], row['QRY_POS'], row['QRY_END'], row['#CHROM'], row['POS'], row['END']
        ))


def check_record_err_string(df, df_qry_fai):
    """
    Runs check_record on each row of `df`, captures exceptions, and returns a Series of error message strings instead
    of failing on the first error. The Series can be added as a column to `df`. For each record where there was no
    error, the field for that record in the returned series is NA (`np.nan`). This function may not be used by the
    pipeline, but is here for troubleshooting alignments.

    :param df: Dataframe of alignment records.
    :param df_qry_fai: Panadas Series with query names as keys and query lengths as values.

    :return: A Series of error messages (or NA) for each record in `df`.
    """

    def _check_record_err_string_check_row(row):
        try:
            check_record(row, df_qry_fai)
            return None
        except Exception as ex:
            return str(ex)

    return df.apply(_check_record_err_string_check_row, axis=1)


def count_cigar(row, allow_m=False):
    """
    Count bases affected by CIGAR operations in an alignment record (row is a Pandas Series from an ailgnment BED).

    Returns a tuple of:
    * ref_bp: Reference bases traversed by CIGAR operations.
    * qry_bp: Query bases traversed by CIGAR operations. Does not include clipped bases.
    * clip_h_l: Hard-clipped bases on the left (upstream) side.
    * clip_s_l: Soft-clipped bases on the left (upstream) side.
    * clip_h_r: Hard-clipped bases on the right (downstream) side.
    * clip_s_r: Soft-clipped bases on the right (downstream) side.

    :param row: Row with CIGAR records as a CIGAR string.
    :param allow_m: If True, allow "M" CIGAR operations. PAV does not allow M operations, this option exists for other
        tools using the PAV library.

    :return: A tuple of (ref_bp, qry_bp, clip_h_l, clip_s_l, clip_h_r, clip_s_r).
    """

    ref_bp = 0
    qry_bp = 0

    clip_s_l = 0
    clip_h_l = 0

    clip_s_r = 0
    clip_h_r = 0

    cigar_list = list(cigar_str_to_tuples(row))

    cigar_n = len(cigar_list)

    index = 0

    while index < cigar_n and cigar_list[index][1] in {'S', 'H'}:
        cigar_len, cigar_op = cigar_list[index]

        if cigar_op == 'S':
            if clip_s_l > 0:
                raise RuntimeError('Duplicate S records (left) at index {}'.format(index))
            clip_s_l = cigar_len

        if cigar_op == 'H':
            if clip_h_l > 0:
                raise RuntimeError('Duplicate H records (left) at index {}'.format(index))

            if clip_s_l > 0:
                raise RuntimeError('S record before H (left) at index {}'.format(index))

            clip_h_l = cigar_len

        index += 1

    while index < cigar_n:
        cigar_len, cigar_op = cigar_list[index]

        # print('{}: {} {}'.format(index, cigar_len, cigar_op))

        if cigar_op in {'=', 'X'}:

            if clip_s_r > 0 or clip_h_r > 0:
                raise RuntimeError(
                    'Found clipped bases before last non-clipped CIGAR operation at operation {} ({}{})'.format(
                        index, cigar_len, cigar_op
                    )
                )

            ref_bp += cigar_len
            qry_bp += cigar_len

        elif cigar_op == 'I':

            if clip_s_r > 0 or clip_h_r > 0:
                raise RuntimeError(
                    'Found clipped bases before last non-clipped CIGAR operation at operation {} ({}{})'.format(
                        index, cigar_len, cigar_op
                    )
                )

            qry_bp += cigar_len

        elif cigar_op == 'D':

            if clip_s_r > 0 or clip_h_r > 0:
                raise RuntimeError(
                    'Found clipped bases before last non-clipped CIGAR operation at operation {} ({}{})'.format(
                        index, cigar_len, cigar_op
                    )
                )

            ref_bp += cigar_len

        elif cigar_op == 'S':

            if clip_s_r > 0:
                raise RuntimeError('Duplicate S records (right) at operation {}'.format(index))

            if clip_h_r > 0:
                raise RuntimeError('H record before S record (right) at operation {}'.format(index))

            clip_s_r = cigar_len

        elif cigar_op == 'H':

            if clip_h_r > 0:
                raise RuntimeError('Duplicate H records (right) at operation {}'.format(index))

            clip_h_r = cigar_len

        elif cigar_op == 'M':

            if not allow_m:
                raise RuntimeError('CIGAR op "M" is not allowed')

            if clip_s_r > 0 or clip_h_r > 0:
                raise RuntimeError(
                    'Found clipped bases before last non-clipped CIGAR operation at operation {} ({}{})'.format(
                        index, cigar_len, cigar_op
                    )
                )

            ref_bp += cigar_len
            qry_bp += cigar_len

        else:
            raise RuntimeError('Bad CIGAR op: ' + cigar_op)

        index += 1

    return ref_bp, qry_bp, clip_h_l, clip_s_l, clip_h_r, clip_s_r


def get_align_bed(align_file, df_qry_fai, hap, min_mapq=0, score_model=None):
    """
    Read alignment file as a BED file that PAV can process. Drops any records marked as unaligned by the SAM flag.

    :param align_file: SAM, CRAM, BAM, anything `pysam.AlignmentFile` can read.
    :param df_qry_fai: Pandas Series with query names as keys and query lengths as values. Index should be cast as
        type str if query names are numeric.
    :param hap: Haplotype assinment for this alignment file (h1 or h2).
    :param min_mapq: Minimum MAPQ. If 0, then all alignments are accepted as long as the unmapped flag is not set.
    :param score_model: Alignment model object (`pavlib.align.score.ScoreModel`) or a configuration string to generate
        a score model object. If `None`, the default score model is used. An alignment score is computed by summing
        the score of each CIGAR operation against this model (match, mismatch, and gap) to create the "SCORE" column.

    :return: BED file of alignment records.
    """

    # Get score model
    score_model = pavlib.align.score.get_score_model(score_model)

    # Get records from SAM
    record_list = list()

    align_index = -1

    with pysam.AlignmentFile(align_file, 'rb') as in_file:
        for record in in_file:

            # Increment align_index
            align_index += 1

            # Skipped unmapped reads
            if record.is_unmapped or record.mapping_quality < min_mapq or len(record.cigar) == 0:
                continue

            # Get length for computing real query positions for rev-complemented records
            qry_len = df_qry_fai[record.query_name]

            # Read tags
            tags = dict(record.get_tags())

            # Determine left hard-clipped bases.
            # pysam query alignment functions are relative to the sequence in the alignment record, not the original
            # sequence. The left-most hard-clipped bases must be added to the query positions to translate to the
            # correct query coordinates (https://github.com/pysam-developers/pysam/issues/1017).
            if len(record.cigartuples) > 0:
                clip_h = record.cigartuples[0][1] if record.cigartuples[0][0] == CIGAR_H else 0
            else:
                clip_h = 0

            cigar_tuples = clip_soft_to_hard(record.cigartuples.copy())

            qry_map_pos = cigar_tuples[0][1] if cigar_tuples[0][0] == CIGAR_H else 0
            qry_map_len = record.query_alignment_end - record.query_alignment_start
            qry_map_end = qry_map_pos + qry_map_len

            if record.query_alignment_start + clip_h != qry_map_pos:
                raise RuntimeError(f'First aligned based from pysam ({record.query_alignment_start}) does not match clipping ({qry_map_pos}) at alignment record {align_index}')

            cigar_string = ''.join(f'{op_len}{CIGAR_CODE_TO_CHAR[op_code]}' for op_code, op_len in cigar_tuples)

            # Disallow alignment match (M) in CIGAR (requires =X for base match/mismatch)
            if CIGAR_M in {op_code for op_code, op_len in cigar_tuples}:
                raise RuntimeError((
                    'Found alignment match CIGAR operation (M) for record {} (Start = {}:{}): '
                    'Alignment requires CIGAR base-level match/mismatch (=X)'
                ).format(record.query_name, record.reference_name, record.reference_start))

            # Save record
            record_list.append(
                pd.Series(
                    [
                        record.reference_name,
                        record.reference_start,
                        record.reference_end,

                        align_index,

                        record.query_name,
                        qry_len - qry_map_end if record.is_reverse else qry_map_pos,
                        qry_len - qry_map_pos if record.is_reverse else qry_map_end,

                        # pavlib.seq.Region(record.query_name, qry_map_pos, qry_map_pos + qry_map_len).to_base1_string(),

                        tags['RG'] if 'RG' in tags else 'NA',
                        tags['AO'] if 'AO' in tags else 'NA',

                        record.mapping_quality,

                        record.is_reverse,
                        f'0x{record.flag:04x}',
                        hap,

                        cigar_string,
                        score_model.score_cigar_tuples(cigar_tuples, rev=True)  # rev for pysam CIGAR (operation code before length)
                    ],
                    index=[
                        '#CHROM', 'POS', 'END',
                        'INDEX',
                        'QRY_ID', 'QRY_POS', 'QRY_END',
                        'RG', 'AO',
                        'MAPQ',
                        'REV', 'FLAGS', 'HAP',
                        'CIGAR', 'SCORE'
                    ]
                )
            )

    # Merge records
    if len(record_list) > 0:
        df = pd.concat(record_list, axis=1).T
    else:
        df = pd.DataFrame(
            [],
            columns=[
                '#CHROM', 'POS', 'END',
                'INDEX',
                'QRY_ID', 'QRY_POS', 'QRY_END',
                'RG', 'AO',
                'MAPQ',
                'REV', 'FLAGS', 'HAP',
                'CIGAR', 'SCORE'
            ]
        )

    # Assign order per query sequence
    df.sort_values(['QRY_ID', 'QRY_POS', 'QRY_END'], inplace=True)

    df['QRY_ORDER'] = -1

    for qry_id in df['QRY_ID'].unique():
        df.loc[df['QRY_ID'] == qry_id, 'QRY_ORDER'] = df.loc[df['QRY_ID'] == qry_id, 'QRY_POS'].rank().astype(int) - 1

    # Reference order
    df.sort_values(['#CHROM', 'POS', 'END', 'QRY_ID'], ascending=[True, True, False, True], inplace=True)

    # Check sanity
    df.apply(check_record, df_qry_fai=df_qry_fai, axis=1)

    # Return BED
    return df


def clip_soft_to_hard(cigar_tuples):
    """
    Strip CIGAR clipping.

    :param cigar_tuples: List of CIGAR operations as tuples (opcode, oplen).

    :return: CIGAR tuples with soft-clipping converted to hard-clipping.
    """

    front_n = 0

    while len(cigar_tuples) > 0 and cigar_tuples[0][0] in {CIGAR_H, CIGAR_S}:
        front_n += cigar_tuples[0][1]
        cigar_tuples = cigar_tuples[1:]

    back_n = 0

    while len(cigar_tuples) > 0 and cigar_tuples[-1][0] in {CIGAR_H, CIGAR_S}:
        back_n += cigar_tuples[-1][1]
        cigar_tuples = cigar_tuples[:-1]

    if len(cigar_tuples) == 0:
        if front_n + back_n == 0:
            raise RuntimeError('Cannot convert soft clipping to hard: No CIGAR records')

        cigar_tuples = [(front_n + back_n, CIGAR_H)]

    else:
        if front_n > 0:
            cigar_tuples = [(CIGAR_H, front_n)] + cigar_tuples

        if back_n > 0:
            cigar_tuples = cigar_tuples + [(CIGAR_H, back_n)]

    return cigar_tuples


def aggregate_alignment_records(df_align, df_qry_fai, score_model=None, min_score=None, colinear_penalty=True, assign_order=True):
    """
    Aggregate colinear alignment records.

    :param df_align: Table of alignment records. MUST be query trimmed (or query- & reference-trimmed)
    :param df_qry_fai: Query FAI.
    :param score_model: Model for scoring INS and DEL between alignment records. If none, use the default model.

    :return:Table of aggregated alignment records.
    """

    # Check parameters
    if score_model is None:
        score_model = pavlib.align.score.get_score_model(pavlib.config.get_config('align_score_model'))

    if min_score is None:
        min_score = score_model.gap(10000)

    min_agg_index = int(10 ** np.ceil(np.log10(
        np.max(df_align['INDEX'])
    ))) - 1  # Start index for aggregated records at the next power of 10

    next_agg_index = min_agg_index + 1

    # Sort
    df_align = df_align.sort_values(['QRY_ID', 'QRY_POS']).copy()
    df_align['INDEX_PREAGG'] = df_align['INDEX']

    # Return existing table if empty
    if df_align.shape[0] == 0:
        df_align['INDEX_PREAGG'] = df_align['INDEX']
        return df_align

    df_align['INDEX_PREAGG'] = df_align['INDEX'].apply(lambda val: [val])
    df_align['MAPQ'] = df_align['MAPQ'].apply(lambda val: [val])
    df_align['FLAGS'] = df_align['FLAGS'].apply(lambda val: [val])

    # Find and aggregate near co-linear records over SVs
    align_records = list()  # Records that were included in a merge

    CLIP_CHAR_SET = {CIGAR_CODE_TO_CHAR[CIGAR_S], CIGAR_CODE_TO_CHAR[CIGAR_H]}

    for qry_id in sorted(set(df_align['QRY_ID'])):
        df = df_align.loc[df_align['QRY_ID'] == qry_id]
        i_max = df.shape[0] - 1

        i = 0
        row1 = df.iloc[i]

        while i < i_max:
            i += 1

            row2 = row1
            row1 = df.iloc[i]

            # Skip if chrom or orientation is not the same
            if row1['#CHROM'] != row2['#CHROM'] or row1['REV'] != row2['REV']:
                align_records.append(row2)
                continue

            # Get reference distance
            if row1['REV']:
                ref_dist = row2['POS'] - row1['END']
            else:
                ref_dist = row1['POS'] - row2['END']

            qry_dist = row1['QRY_POS'] - row2['QRY_END']

            if qry_dist < 0:
                raise RuntimeError(f'Query distance is negative: {qry_dist}: alignment indexes {row1["INDEX"]} and {row2["INDEX"]}')

            if ref_dist >= 0:
                # Contiguous in reference space, check query space

                # Score gap between the alignment records
                score = score_model.gap(ref_dist) + score_model.gap(qry_dist)

                score_gap = score + (
                    score_model.gap(np.abs(qry_dist - ref_dist)) if colinear_penalty else 0
                )

                if score_gap < min_score:
                    align_records.append(row2)
                    continue

                #
                # Aggregate
                #
                row1 = row1.copy()

                # Set query position
                row1['QRY_POS'] = row2['QRY_POS']

                # Get rows in order
                if row1['REV']:
                    row_l = row1
                    row_r = row2

                    row1['END'] = row2['END']

                    row1['TRIM_REF_R'] = row2['TRIM_REF_R']
                    row1['TRIM_QRY_R'] = row2['TRIM_QRY_R']

                else:
                    row_l = row2
                    row_r = row1

                    row1['POS'] = row2['POS']

                    row1['TRIM_REF_L'] = row2['TRIM_REF_L']
                    row1['TRIM_QRY_L'] = row2['TRIM_QRY_L']

                # Set records
                row1['FLAGS'] = row1['FLAGS'] + row2['FLAGS']
                row1['MAPQ'] = row1['MAPQ'] + row2['MAPQ']
                row1['INDEX_PREAGG'] = row1['INDEX_PREAGG'] + row2['INDEX_PREAGG']
                row1['SCORE'] = row1['SCORE'] + row2['SCORE']

                if 'RG' in row1:
                    row1['RG'] = np.nan

                if 'AO' in row1:
                    row1['AO'] = np.nan

                # Merge CIGAR strings
                cigar_l = list(cigar_str_to_tuples(row_l['CIGAR']))
                cigar_r = list(cigar_str_to_tuples(row_r['CIGAR']))

                while len(cigar_l) > 0 and cigar_l[-1][1] in CLIP_CHAR_SET:  # Tail of left record
                    cigar_l = cigar_l[:-1]

                while len(cigar_r) > 0 and cigar_r[0][1] in CLIP_CHAR_SET:  # Head of right record
                    cigar_r = cigar_r[1:]

                ins_len = qry_dist
                del_len = ref_dist

                if qry_dist > 0:
                    while len(cigar_l) > 0 and cigar_l[-1][1] == CIGAR_CODE_TO_CHAR[CIGAR_I]:  # Concat insertions (no "...xIxI..." in CIGAR)
                        ins_len += cigar_l[-1][0]
                        cigar_l = cigar_l[:-1]

                    cigar_l.append((ins_len, CIGAR_CODE_TO_CHAR[CIGAR_I]))

                if ref_dist > 0:
                    while len(cigar_l) > 0 and cigar_l[-1][1] == CIGAR_CODE_TO_CHAR[CIGAR_D]:  # Concat deletions (no "...xDxD..." in CIGAR)
                        del_len += cigar_l[-1][0]
                        cigar_l = cigar_l[:-1]

                    cigar_l.append((del_len, CIGAR_CODE_TO_CHAR[CIGAR_D]))

                row1['CIGAR'] = ''.join([f'{record[0]}{record[1]}' for record in cigar_l + cigar_r])

                row1['SCORE'] = score_model.score_cigar_tuples(row1['CIGAR'])

                # Set alignment indexes
                if row2['INDEX'] < min_agg_index:
                    # Use the next aggregate index
                    row1['INDEX'] = next_agg_index
                    next_agg_index += 1

                else:
                    # row2 was aggregated, use its aggregate index
                    row1['INDEX'] = row2['INDEX']

                # Check
                check_record(row1, df_qry_fai)

            else:
                align_records.append(row2)

        # Add last record
        align_records.append(row1)

    # Concatenate records
    df = pd.concat(align_records, axis=1).T

    df['MAPQ'] = df['MAPQ'].apply(lambda val: ','.join([str(v) for v in val]))
    df['FLAGS'] = df['FLAGS'].apply(lambda val: ','.join([str(v) for v in val]))
    df['INDEX_PREAGG'] = df['INDEX_PREAGG'].apply(lambda val: ','.join([str(v) for v in val]))

    # Assign order per query sequence
    if assign_order:
        df.sort_values(['QRY_ID', 'QRY_POS', 'QRY_END'], inplace=True)

        df['QRY_ORDER'] = -1

        for qry_id in df['QRY_ID'].unique():
            df.loc[df['QRY_ID'] == qry_id, 'QRY_ORDER'] = df.loc[df['QRY_ID'] == qry_id, 'QRY_POS'].rank().astype(int) - 1

    # Reference order
    df.sort_values(['#CHROM', 'POS', 'END', 'QRY_ID'], ascending=[True, True, False, True], inplace=True)

    # Check sanity (whole table, modified records already checked, should pass)
    df.apply(check_record, df_qry_fai=df_qry_fai, axis=1)

    return df
