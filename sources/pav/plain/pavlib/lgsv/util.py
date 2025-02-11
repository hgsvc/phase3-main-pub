"""
Variant call utilities.
"""

import numpy as np
import pandas as pd

import Bio.Seq
import pavlib
import svpoplib


class CallerResources(object):
    """
    Container of resources needed by routines attempting to resolve variantns.

    Attributes:
        * df_align_qry: Alignment table (QRY trimmed)
        * df_align_qryref: Alignment table (QRY & REF trimmed)
        * df_align_none: Alignment table (No trimming)
        * qry_fa_name: Query FASTA filename.
        * ref_fa_name: Reference FASTA filename.
        * hap: Haplotype name.
        * score_model: Alignemnt score model.
        * k_util: K-mer utility.
        * inv_params: Inversion parameters.
        * kde: KDE model.

        * align_lift: Object for lifting alignment coordinates between query and reference through the alignment.
        * cache_qry_upper: Query sequence cache used for homology searches. Caches the last query sequence in
            upper-case.
        * cache_ref_upper: Reference sequence cache used for homology searches. Caches the last reference sequence in
            upper-case.
        * off_gap_mult: Off-target gap multiplier (e.g. reference gap for an insertion or query gap for a deletion).
        *
    """

    def __init__(self,
                 df_align_qry,
                 df_align_qryref,
                 df_align_none,
                 qry_fa_name, ref_fa_name,
                 hap,
                 score_model=None, k_util=None,
                 inv_params=None,
                 kde=None,
                 log_file=None,
                 config_params=None
                 ):

        if config_params is None:
            config_params = pavlib.config.get_config_dict()

        self.df_align_qry = df_align_qry
        self.df_align_qryref = df_align_qryref
        self.df_align_none = df_align_none

        self.qry_fa_name = qry_fa_name
        self.qry_fai_name = self.qry_fa_name + '.fai'

        self.ref_fa_name = ref_fa_name
        self.ref_fai_name = self.ref_fa_name + '.fai'

        self.ref_fai = svpoplib.ref.get_df_fai(self.ref_fai_name)
        self.qry_fai = svpoplib.ref.get_df_fai(self.qry_fai_name)

        self.hap = hap

        if k_util is None:
            k_util = pavlib.kmer.KmerUtil(config_params.inv_k_size)

        self.k_util = k_util

        if score_model is None:
            score_model = pavlib.align.score.get_score_model()

        self.score_model = score_model

        self.align_lift = pavlib.align.lift.AlignLift(self.df_align_qry, self.qry_fai)

        self.cache_qry_upper = pavlib.lgsv.util.SeqCache(self.qry_fa_name, uppercase=True)  # uppercase for homology searches
        self.cache_ref_upper = pavlib.lgsv.util.SeqCache(self.ref_fa_name, uppercase=True)

        if kde is None:
            kde = pavlib.kde.KdeTruncNorm()

        self.kde = kde

        if inv_params is None:
            inv_params = dict()
        else:
            inv_params = dict(inv_params)

        for params in list(inv_params.keys()):
            if params not in {'nc_ref', 'nc_qry', 'region_limit', 'init_expand', 'min_kmers',
                              'max_ref_kmer_count', 'repeat_match_prop', 'min_inv_kmer_run', 'min_qry_ref_prop'}:
                del inv_params[params]

        self.inv_params = inv_params

        self.log_file = log_file
        self.config_params = config_params


def record_to_paf(row_seg, ref_fai, qry_fai, mapq_summary='max'):
    """
    Convert the row of a segment table to PAF format.

    `row_seg` is a complex segment record with "MAPQ" and "CIGAR" fields added.

    :param row_seg: Segment table row.
    :param ref_fai: Reference FASTA index.
    :param qry_fai: Query FASTA index.
    :param mapq_summary: If multiple alignment records were aggregated, then the MAPQ value is a list of MAPQ values
        from the original alignments. When multiple MAPQ vaules are found, summarize them to a single value with this
        approach. "max": maximum value (default), "min": minimum value, "mean": average value.

    :return: PAF record row.
    """

    match_n = 0
    align_len = 0
    cigar_index = -1

    cigar_list = list(pavlib.align.util.cigar_str_to_tuples(row_seg['CIGAR']))

    # Remove clipping and adjust coordinates
    if cigar_list[0][1] == 'H':
        cigar_list = cigar_list[1:]
        cigar_index += 1

    if cigar_list[0][1] == 'S':
        cigar_list = cigar_list[1:]
        cigar_index += 1

    if cigar_list[-1][1] == 'H':
        cigar_list = cigar_list[:-1]

    if cigar_list[-1][1] == 'S':
        cigar_list = cigar_list[:-1]

    cigar = ''.join([f'{op_len}{op_code}' for op_len, op_code in cigar_list])

    # Process CIGAR operations
    for op_len, op_code in cigar_list:
        cigar_index += 1

        if op_code == '=':
            match_n += op_len
            align_len += op_len

        elif op_code in {'X', 'I', 'D'}:
            align_len += op_len

        elif op_code in {'H', 'S'}:
            raise RuntimeError(f'Unhandled clipping in CIGAR string: {op_code} at CIGAR index {cigar_index}: Expected clipped bases at the beginning and end of the CIGAR string only.')

        else:
            raise RuntimeError(f'Unrecognized CIGAR op code: {op_code} at CIGAR index {cigar_index}')

    # Set strand
    if 'STRAND' in row_seg:
        strand = row_seg['STRAND']
    elif 'REV' in row_seg:
        strand = '-' if row_seg['REV'] else '+'
    elif 'IS_REV' in row_seg:
        strand = '-' if row_seg['IS_REV'] else '+'
    else:
        raise RuntimeError(f'Missing "STRAND", "REV", or "IS_REV" column in segment table: Record {row_seg["INDEX"] if "INDEX" in row_seg else row_seg.name}')

    # Adjust MAPQ (might be a list of MAPQ values)
    if isinstance(row_seg['MAPQ'], str):
        mapq_list = [int(v) for v in row_seg['MAPQ'].split(',')]

        if mapq_summary == 'max':
            mapq = np.max(mapq_list)
        elif mapq_summary == 'min':
            mapq = np.min(mapq_list)
        elif mapq_summary == 'mean':
            mapq = np.mean(mapq_list)
        else:
            raise RuntimeError(f'Unrecognized mapq_summary: {mapq_summary}')
    else:
        mapq = row_seg['MAPQ']

    # Create PAF record
    return pd.Series(
        [
            row_seg['QRY_ID'],
            qry_fai[row_seg['QRY_ID']],
            row_seg['QRY_POS'],
            row_seg['QRY_END'],
            strand,
            row_seg['#CHROM'],
            ref_fai[row_seg['#CHROM']],
            row_seg['POS'],
            row_seg['END'],
            match_n,
            align_len,
            mapq,
            cigar
        ],
        index=[
            'QRY_NAME',
            'QRY_LEN',
            'QRY_POS',
            'QRY_END',
            'STRAND',
            'CHROM',
            'CHROM_LEN',
            'CHROM_POS',
            'CHROM_END',
            'MISMATCH_N',
            'ALIGN_BLK_LEN',
            'MAPQ',
            'CIGAR'
        ]
    )


class SeqCache:
    """
    Keep a cache of a sequence string in upper-case. Stores the last instance of the sequence and the ID. When a
    new ID is requested, the old sequnece is discarded and the new one is loaded.
    """

    def __init__(self, fa_filename, uppercase=True):
        """
        Create a cache object to read from indexed FASTA file `fa_filename`.

        :param fa_filename: Indexed FASTA file name.
        :param uppercase: `True` if sequences should be made upper-case, otherwise, preserve case.
        """

        self.fa_filename = fa_filename
        self.uppercase = uppercase

        self.id = None
        self.seq = None
        self.is_rev = None

    def get(self, sequence_id, is_rev):
        """
        Get a sequence. Returns the cached version if ID matches the current ID, otherwise, the correct sequence is
        retrieved, cached, and returned.

        :param sequence_id: Sequence ID string or Region.
        :param is_rev: `True` if the sequence is reverse-complemented. Retrieving the same sequence ID as the cached
            sequence with `is_rev` mismatch will reload the sequence in the requested orientation.

        :return: Sequence.
        """

        if self.id != sequence_id or self.is_rev != is_rev:
            new_seq = pavlib.seq.region_seq_fasta(
                sequence_id, self.fa_filename
            )

            if self.uppercase:
                new_seq = new_seq.upper()

            if is_rev:
                new_seq = str(Bio.Seq.Seq(new_seq).reverse_complement())

            self.seq = new_seq

            self.id = sequence_id
            self.is_rev = is_rev

        return self.seq


def dot_graph_writer(out_file, df_align, chain_set, optimal_interval_list, sv_dict, graph_name='Unnamed_Graph', forcelabels=True, anchor_width=2.5):

    # Header
    out_file.write(f'graph {graph_name} {{\n')

    # Attributes
    if forcelabels:
        out_file.write('    forcelabels=true;\n')

    out_file.write('    overlap=false;\n')

    # Anchor and interval sets
    optimal_interval_set = set(optimal_interval_list)
    variant_interval_set = set(sv_dict.keys())

    optimal_anchor_set = {index for index_pair in optimal_interval_set for index in index_pair}
    variant_anchor_set = {index for index_pair in sv_dict.keys() for index in index_pair}

    # optimal_interval_set = set(optimal_interval_list)

    # Add nodes
    for index, row in df_align.iterrows():

        if index in optimal_anchor_set:
            color = 'blue'
        elif index in variant_anchor_set:
            color = 'black'
        else:
            color = 'gray33'

        width = anchor_width if index in variant_anchor_set else 1

        out_file.write(f'    n{index} [label="{index} ({row["INDEX"]}) - {row["#CHROM"]}:{row["POS"]}-{row["END"]} {"-" if row["REV"] else "+"} s={row["SCORE"]}", penwidth={width}, color="{color}"]\n')

    # Add candidate edges
    for start_index, end_index in chain_set:

        if (start_index, end_index) in optimal_interval_set:
            color = 'blue'
        elif (start_index, end_index) in variant_interval_set:
            color = 'black'
        else:
            color = 'gray33'

        width = anchor_width if (start_index, end_index) in variant_interval_set else 1

        if (start_index, end_index) not in sv_dict or sv_dict[start_index, end_index].is_null():
            var_name = 'NullVar'
        elif sv_dict[start_index, end_index].is_patch:
            var_name = 'AlignPatch (No Variant)'
        else:
            var_name = sv_dict[start_index, end_index].variant_id

        out_file.write(f'    n{start_index} -- n{end_index} [label="{var_name} (s={sv_dict[start_index, end_index].score_variant})", penwidth={width}, color="{color}"]\n')

    # Add adjacent edges (not anchor candidates)
    for start_index in range(df_align.shape[0] - 1):
        end_index = start_index + 1

        if (start_index, end_index) in optimal_interval_set:
            color = 'blue'
        elif (start_index, end_index) in variant_interval_set:
            color = 'black'
        else:
            color = 'gray33'

        if (start_index, end_index) not in chain_set:
            out_file.write(f'    n{start_index} -- n{end_index} [penwidth=1, color="{color}"]\n')

    # Done
    out_file.write('}\n')


def get_min_anchor_score(min_anchor_score, score_model):

    if isinstance(min_anchor_score, str):
        min_anchor_score_str = min_anchor_score.strip()

        if len(min_anchor_score_str) == 0:
            raise ValueError('min_anchor_score parameter is empty')

        if min_anchor_score_str.lower().endswith('bp'):
            try:
                min_anchor_score_bp = abs(float(min_anchor_score_str[:-2].strip()))
            except ValueError:
                raise ValueError(f'min_anchor_score is a number before "bp": {min_anchor_score}')

            if min_anchor_score_bp <= 0.0:
                return 0.0

            return score_model.match(min_anchor_score_bp)

        else:
            try:
                return abs(float(min_anchor_score))
            except ValueError:
                raise ValueError(f'min_anchor_score is a string that does not represent a numeric value: {min_anchor_score}')

    else:
        try:
            # noinspection PyTypeChecker
            return float(min_anchor_score)
        except ValueError:
            raise ValueError(f'min_anchor_score is not a string or numeric: type={type(min_anchor_score)}')
