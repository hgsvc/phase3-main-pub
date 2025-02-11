"""
Large variant discovery call methods
"""

import abc
import numpy as np
import pandas as pd
import scipy.stats

import pavlib

#
# Variant call objects
#

CALL_SOURCE = 'ALNTRUNC'

REF_TRACE_COLUMNS_PRE = ['#CHROM', 'POS', 'END', 'DEPTH', 'QRY_ID', 'INDEX']
REF_TRACE_COLUMNS_POST = ['FWD_COUNT', 'REV_COUNT', 'TYPE', 'LEN']
REF_TRACE_COLUMNS = REF_TRACE_COLUMNS_PRE + REF_TRACE_COLUMNS_POST


#
# Variant discovery routines
#

def try_variant(var_type, interval, caller_resources, best_variant, var_region_kde):
    """
    Try calling a variant of type `var_type` (a subclass of `Variant`). Check the variant against `best_variant` and
    return the variant of these two (new variant and best variant) with the highest score.

    :param var_type: Variant call type to try (A subclass of `Variant`).
    :param interval: Alignment interval.
    :param caller_resources: Caller resources.
    :param best_variant: Best variant call so far. If there is no variant to check against, a `NullVariant` object
        should be used as this parameter.
    :param var_region_kde: For regions where the inserted and deleted segments are of similar size, this object
        describes how well the inserted sequence matches the deleted sequence in forward or reverse orientation. If
        the forward match is too high, then do not attempt a variant call.

    :return: Best variant call between the new variant of type `var_type` and `best_variant`. May return a `NullVariant`
        if both `best_variant` is a `NullVariant` and the interval does not match a variant of type `var_type`.
    """

    if best_variant is None:
        raise ValueError('best_variant cannot be None')

    # Try variant call
    try:
        variant = var_type(interval, caller_resources, var_region_kde)

    except VariantMatchError as e:
        # Null variant is this locus does not match a variant of this type
        variant = NullVariant()

    if variant.score_variant <= best_variant.score_variant:
        return best_variant

    return variant

def score_segment_transitions(df_segment, caller_resources):
    """
    Score complex variant on template switches and gaps for each segment.

    :param df_segment: Segment table.
    :param caller_resources: Caller resources.

    :return: Variant score.
    """

    score_variant = caller_resources.score_model.template_switch() * \
                    (df_segment.shape[0] - 1)  # Template switches between segments (n + 2 alignment records (n segments + 2 anchors), n - 1 template switches (between each segment including each anchor)

    for i in range(1, df_segment.shape[0] - 1):  # Gap per segment
        score_variant += caller_resources.score_model.gap(
            df_segment.loc[i, 'LEN_QRY']
        )

    return score_variant


#
# Variant call objects
#

class Variant(object, metaclass=abc.ABCMeta):
    """
    Base class for variant call objects.

    Attributes:
        vartype: Variant class
        svtype: Variant type
        svsubtype: Variant subtype (e.g. TANDEM)
        svlen: Variant length
        call_source: Type of evidence supporting the variant call
        seq: Variant sequence or "*" if the sequence is not set.
        score_variant: Variant score

        hom_ref: Breakpoint homology in the reference (set by complete_anno())
        hom_qry: Breakpoint homology in the query (set by complete_anno())

        interval: Variant interval
        df_ref_trace: Variant reference trace

        is_complete_anno: Set when variant annotations are complete (only on accepted variants)

    Attributes derived on demand or through the interval:
        variant_id: Variant ID
        chrom: Reference region chrom
        pos: Reference region position
        end: Reference region end
        qry_id: Query ID
        qry_pos: Query region position
        qry_end: Query region end
        is_rev: Is rev
        strand: Strang
        chain_start_index: Start index from the alignment chain
        chain_end_index: End index from the alignment chain
        region_ref: Reference region from the interval
        region_qry: Query region from the interval
    """

    def __init__(self, interval):

        # Set null type
        self.interval = interval

        self.score_variant = -np.inf

        self.df_ref_trace = None

        self.vartype = 'NullVariant'
        self.svtype = 'NullVariant'
        self.svsubtype = None
        self.svlen = 0

        self.hom_ref = None
        self.hom_qry = None

        self.is_complete_anno = False

        self.call_source = CALL_SOURCE
        self.seq = '*'

        self.is_complete_anno = False

        self.region_qry = None

        # Check for null type
        if interval is None:

            self.chrom = '*'
            self.pos = 0
            self.end = 0

            self.qry_id = '*'
            self.qry_pos = None
            self.qry_end = None

            self.anchor_score_min = -np.inf
            self.anchor_score_max = np.inf

            self.is_complete_anno = True

            return

        self.chrom = interval.region_ref.chrom
        self.pos = None
        self.end = None

        self.qry_id = interval.region_qry.chrom
        self.qry_pos = None
        self.qry_end = None

        self.anchor_score_min = np.min(self.interval.df_segment.loc[self.interval.df_segment['IS_ANCHOR'], 'SCORE'])
        self.anchor_score_max = np.max(self.interval.df_segment.loc[self.interval.df_segment['IS_ANCHOR'], 'SCORE'])

        self.df_kde = None

        self.is_patch = False  # Patch variant, represents a bridge across alignment artifacts (no real variant)

    def __getattr__(self, name):
        """
        Get additional attributes by name.

        :param name: Attribute name.

        :return: Attribute value.
        """

        if name == 'variant_id':
            return f'{self.chrom}-{self.pos}-{str(self.svtype).upper()}-{self.svlen}'

        if name == 'is_rev':
            return self.interval.is_rev if not self.is_null() else False

        elif name == 'strand':
            return ('-' if self.is_rev else '+') if not self.is_null() else '*'

        elif name == 'chain_start_index':
            return self.interval.chain_node.start_index if not self.is_null() else '*'

        elif name == 'chain_end_index':
            return self.interval.chain_node.end_index if not self.is_null() else '*'

        elif name == 'region_ref':
            return self.interval.region_ref if not self.is_null() else pavlib.seq.Region('*', 0, 1)

        elif name == 'QRY_REGION':
            return self.interval.region_qry if not self.is_null() else pavlib.seq.Region('*', 0, 1)

        elif name == 'seg_n':
            return self.interval.df_segment.shape[0] if not self.is_null() else 0

        raise AttributeError(f'Variant has no attribute: {name}')

    def __repr__(self):
        """
        String representation of a variant object.

        :return: String representation.
        """

        return f'Variant({self.variant_id}, ref={self.region_ref}, qry={self.region_qry}, strand={self.strand}, score={self.score_variant})'

    def complete_anno(self):
        """
        Complete annotations on this variant call.

        Initial variant calls are not completed so they can be prioritized
        before spending the CPU cycles and IO time pulling sequences to complete annotations.

        Implementations of this method should complete annotations and set `self.is_complete_anno` to `True`.
        """

        if self.is_complete_anno or self.is_null():
            return

        # Variant type implementation
        self.complete_anno_impl()
        self.is_complete_anno = True

    def is_null(self):
        """
        Determine if variant is a null call.

        :return: `True` if the variant call is Null (does not represent a real variant), and 'False' if the variant call
            exists.
        """

        return self.interval is None

    def complete_anno_impl(self):
        """
        Complete annotations. To be implemented by the variant subclass.
        """
        pass

    def row(self):
        """
        Get a Pandas Series object representing a variant call table row for this variant call.

        Must be called after `complete_anno()`.

        :return: Pandas Series object of this variant call.
        """

        self.complete_anno()

        return self.row_impl()

    @abc.abstractmethod
    def row_impl(self):
        """
        Variant type implementation of row().

        :return: A Pandas Series object representing a variant call table row for this variant call.
        """
        raise NotImplementedError


# TODO: correct untemplated insertion breakpoints for homology at breakpoints
#
# Find untemplated insertion
#
# Trim reference if ends overlap around an untemplated insertion. This is likely homology at the breakpoint
# where the homologous sequence at each end of the insertion was aligned to the reference.
#
#                             [    HOM REGION    ]
#     Ref: --------->--------->--------->--------->--------->--------->--------->
# Align 1: --------->--------->--------->--------->
# Align 2:                     --------->--------->--------->--------->--------->
#
# Alignments 1 & 2 skip query sequence between them (untemplated insertion). This does not occur for strict templated
# insertions because query trimming will remove redundancy.
class InsertionVariant(Variant):
    """
    Insertion variant call.
    """

    # Simple insertion (INS = unaligned or aligned to another site)
    # INS:                           -------------------
    #                                        ||
    # Qry: --->--------->---------->---------> --------->-------->--------->--------->-----
    # Ref: --------------------------------------------------------------------------------

    def __init__(self, interval, caller_resources, var_region_kde=None):
        Variant.__init__(self, interval)

        # Null variant (for creating table headers when no variants are found)
        if interval is None:
            return

        # Return immediately to leave the variant call a Null type
        if interval.seg_n != 1 or interval.len_qry == 0:
            return

        # Reference gap penalty: Penalize insertions for overlaps or deletions in the reference at the INS breakpoint
        len_ref = np.abs(self.interval.len_ref)

        ref_overlap = len_ref if self.interval.len_ref < 0 else 0
        self.svlen = len(interval.region_qry) + ref_overlap

        self.score_variant = \
            caller_resources.score_model.gap(self.svlen) + \
            caller_resources.score_model.gap(abs(len_ref)) * \
                caller_resources.config_params.lg_off_gap_mult

        self.pos = self.region_ref.pos
        self.end = self.pos + 1
        self.vartype = 'SV' if self.svlen > 50 else 'INDEL'
        self.svtype = 'INS'

        self.region_qry = interval.region_qry

    def row_impl(self):
        return pd.Series(
            [
                self.chrom, self.pos, self.end, self.variant_id, self.svtype, self.svlen,
                str(self.region_qry), self.strand,
                self.svsubtype, self.score_variant,
                self.call_source,
                self.anchor_score_min, self.anchor_score_max,
                f'{self.chain_start_index}-{self.chain_end_index}'
            ],
            index=[
                '#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN',
                'QRY_REGION', 'QRY_STRAND',
                'SVSUBTYPE', 'VAR_SCORE',
                'CALL_SOURCE',
                'ANCHOR_SCORE_MIN', 'ANCHOR_SCORE_MAX',
                'INTERVAL'
            ]
        )


class TandemDuplicationVariant(InsertionVariant):
    """
    Tandem duplication variant call.
    """

    # Tandem duplication (TD)
    # Repeats:                                [> REP >]            [> REP >]
    # Qry 1:    --------->--------->--------->--------->--------->--------->
    # Qry 2:                                  --------->--------->--------->--------->--------->--------->
    #
    # Directly-oriented repeats may mediated tandem repeats. Look at alignment-trimming in three ways:
    # * Qry trimming: Identifies TD if redundant query bases are removed, but queries still overlap
    # * Qry & Ref trimming: Find a breakpoint for an insertion call.
    # * No trimming: Identify the repeats at the ends of the TD.
    #
    # The resulting call is an INS call with a breakpoint placed in a similar location if the TD was called as an
    # insertion event in the CIGAR string. The variant is annotated with the DUP locus, and the sequence of both copies
    # is

    def __init__(self, interval, caller_resources, var_region_kde=None):
        Variant.__init__(self, interval)

        if interval is None:
            return

        if interval.seg_n != 0 or interval.len_ref >= 0:
            return

        if interval.df_segment.iloc[0]['INDEX'] not in caller_resources.df_align_qryref:
            return

        if interval.df_segment.iloc[-1]['INDEX'] not in caller_resources.df_align_qryref:
            return

        # Determine left-most breakpoint using alignment trimming for homology
        self.pos = caller_resources.df_align_qryref.loc[interval.df_segment.iloc[0]['INDEX'], 'END']
        self.end = self.pos + 1

        qry_pos = caller_resources.df_align_qryref.loc[interval.df_segment.iloc[0]['INDEX'], 'QRY_END']
        qry_end = caller_resources.df_align_qryref.loc[interval.df_segment.iloc[-1]['INDEX'], 'QRY_POS']

        self.svlen = qry_end - qry_pos

        self.score_variant = \
            caller_resources.score_model.gap(self.svlen)

        self.vartype = 'SV' if self.svlen > 50 else 'INDEL'
        self.svtype = 'INS'
        self.svsubtype = 'TD'

        self.region_qry = interval.region_qry


class DeletionVariant(Variant):
    """
    Deletion variant call.
    """

    # Simple deletion
    # Qry: ---->---------->--------->                  --------->-------->--------->-------
    # Ref: --------------------------------------------------------------------------------

    def __init__(self, interval, caller_resources, var_region_kde=None):
        Variant.__init__(self, interval)

        # Null variant (for creating table headers when no variants are found)
        if interval is None:
            return

        # Return immediately to leave the variant call a Null type
        if interval.len_ref <= 0:
            return

        self.pos = self.interval.region_ref.pos
        self.end = self.interval.region_ref.end
        self.svlen = len(self.interval.region_ref)

        self.region_qry = interval.region_qry

        self.score_variant = caller_resources.score_model.gap(self.svlen) + \
            caller_resources.score_model.gap(abs(self.interval.len_qry)) * \
              caller_resources.config_params.lg_off_gap_mult

        self.vartype = 'SV' if self.svlen > 50 else 'INDEL'
        self.svtype = 'DEL'

    def row_impl(self):
        return pd.Series(
            [
                self.chrom, self.pos, self.end, self.variant_id, self.svtype, self.svlen,
                str(self.region_qry), self.strand,
                self.score_variant,
                self.call_source,
                self.anchor_score_min, self.anchor_score_max,
                f'{self.chain_start_index}-{self.chain_end_index}'
            ],
            index=[
                '#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN',
                'QRY_REGION', 'QRY_STRAND',
                'VAR_SCORE',
                'CALL_SOURCE',
                'ANCHOR_SCORE_MIN', 'ANCHOR_SCORE_MAX',
                'INTERVAL'
            ]
        )


class InversionVariant(Variant):
    """
    Balanced inversion variant call.
    """

    def __init__(self, interval, caller_resources, var_region_kde=None):
        Variant.__init__(self, interval)

        # Note: Breakpoints may not be placed consistently in inverted repeats on each end, the aligner makes an
        # alignment decision independently for each breakpoint. Therefore, the inverted segment may not align to the
        # reference gap between the repeats. An alignment penalty should be applied to discourage inversion calls
        # that do not fit the classic model (an inverted alignment in a reference gap), but not penalize alignment
        # this common alignment artifact.
        #
        #            [ > >  Repeat  > > ]                                        [ < <  Repeat  < < ]
        # Flank: --->--------->------->                                                          --->----
        #   INV:       <---------<---------<---------<---------<---------<---------
        # Trace:  NML    |   DUP      |                  INV                       |    DEL      |   NML
        #
        # The pattern above is apparent in trimmed alignments (by query trimming). The un-trimmed alignment will
        # typically align through both repeats:
        #
        #            [ > >  Repeat  > > ]                                        [ < <  Repeat  < < ]
        # Flank: --->--------->------->--                                        --------->------->------
        #   INV:     --<---------<---------<---------<---------<---------<---------<---------<-------
        # Trace:  NML    |   DUP      |                  INV                       |    DEL      |   NML
        #
        # These two sources of information are combined. The qry-trimmed alignment is used to score the inversion after
        # redundant alignments are removed so complex patterns can be penalized appropriately. The untrimmed alignment
        # is used for setting inner- and outer-breakpoint locations.

        self.region_ref_inner = None  # Non-None is a sentinel for a found variant
        self.region_qry_inner = None
        self.region_ref_outer = None
        self.region_qry_outer = None

        self.call_subtype = 'NA'
        self.size_gap = 0

        # Base inversion checks
        if interval is None or interval.len_ref <= 0 or interval.len_qry <= 0:
            # Balanced inversions consume ref & query bases

            return

        if min(interval.len_ref, interval.len_qry) / max(interval.len_ref, interval.len_qry) < 0.5:
            # Neither ref or query region may be 2x greater than the other

            return

        df_int = interval.df_segment.loc[
            ~ interval.df_segment['IS_ANCHOR']
        ]

        is_prox = (
            df_int['IS_ALIGNED'] &
            (df_int['POS'] < interval.region_ref.end) &
            (df_int['END'] > interval.region_ref.pos)
        )

        is_dist = (~ is_prox) & df_int['IS_ALIGNED']

        if np.sum(df_int.loc[is_dist, 'LEN_QRY']) > interval.len_qry * 0.1:
            # No more than 10% aligned outside of the inversion site (allow for small alignments)
            return

        if not var_region_kde.kde_inv and np.sum(df_int.loc[is_prox & (df_int['IS_REV'] == interval.is_rev), 'LEN_QRY']) > np.sum(df_int.loc[is_prox & (df_int['IS_REV'] != interval.is_rev), 'LEN_QRY']):
            # Inveted by KDE, or Inverted segment lengths should outweigh the non-inverted segments if the inversion is called by alignment.

            return

        # Call variant by KDE
        if var_region_kde is not None and var_region_kde.kde_inv:

            # Set region boundaries
            self.region_ref_inner = interval.region_ref
            self.region_qry_inner = interval.region_qry

            self.region_ref_outer = self.region_ref_inner
            self.region_qry_outer = self.region_qry_inner

            self.call_subtype = 'KDE'

        # Try by alignment
        local_inv = False

        if self.region_ref_inner is None:
            inv_ref_flank = len(interval.region_qry) * 0.1  # Allow inverted segment to be within 5% of the reference region

            df_int = interval.df_segment.loc[
                ~ interval.df_segment['IS_ANCHOR']
            ]

            if df_int.shape[0] == 1:
                row_inv = df_int.iloc[0]

                local_inv = (
                    row_inv['#CHROM'] == interval.chrom
                ) and (
                    row_inv['POS'] < interval.region_ref.pos + inv_ref_flank
                ) and (
                    row_inv['END'] > interval.region_ref.end - inv_ref_flank
                ) and (
                    row_inv['IS_REV'] != interval.is_rev
                )

        if local_inv:
            self.region_ref_inner = interval.region_ref
            self.region_qry_inner = interval.region_qry

            self.region_ref_outer = self.region_ref_inner
            self.region_qry_outer = self.region_qry_inner

            self.call_subtype = 'ALN'

        # Call variant
        if self.region_ref_inner is not None:

            self.size_gap = np.abs(len(self.region_ref_inner) - len(self.region_qry_inner))
            self.svlen = len(self.region_ref_inner)

            self.score_variant = (
                caller_resources.score_model.gap(self.size_gap) +  # Penalize differences between reference gap and inverted alignment
                caller_resources.score_model.template_switch() * 2
            )

            self.chrom = self.region_ref_inner.chrom
            self.pos = self.region_ref_inner.pos
            self.end = self.region_ref_inner.end

            self.region_qry = self.region_qry_inner

            self.vartype = 'SV'
            self.svtype = 'INV'

        return

    def row_impl(self):

        return pd.Series(
            [
                self.chrom, self.pos, self.end, self.variant_id, self.svtype, self.svlen,
                str(self.region_qry_inner), self.strand,
                self.score_variant,
                str(self.region_ref_outer), str(self.region_qry_outer),
                self.size_gap,
                f'{self.call_source}:{self.call_subtype}',
                self.anchor_score_min, self.anchor_score_max,
                f'{self.chain_start_index}-{self.chain_end_index}'
            ],
            index=[
                '#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN',
                'QRY_REGION', 'QRY_STRAND',
                'VAR_SCORE',
                'REGION_REF_OUTER', 'REGION_QRY_OUTER',
                'ALIGN_SIZE_GAP',
                'CALL_SOURCE',
                'ANCHOR_SCORE_MIN', 'ANCHOR_SCORE_MAX',
                'INTERVAL'
            ]
        )


class ComplexVariant(Variant):
    """
    Complex variant call.
    """

    def __init__(self, interval, caller_resources, var_region_kde=None):
        Variant.__init__(self, interval)

        # Null variant (for creating table headers when no variants are found)
        if interval is None or interval.len_qry <= 0:
            self.struct_qry = None
            self.struct_ref = None

            return

        # Get reference trace
        self.df_ref_trace = get_reference_trace(self.interval)

        # Compute variant score
        self.score_variant = \
            score_segment_transitions(interval.df_segment, caller_resources)

        for svlen in self.df_ref_trace.loc[self.df_ref_trace['TYPE'] == 'DEL', 'LEN']:
            self.score_variant += caller_resources.score_model.gap(svlen)

        self.vartype = 'SV'
        self.svtype = 'CPX'
        self.svlen = len(interval.region_qry) + np.abs(interval.len_ref)

        self.struct_qry = None
        self.struct_ref = None
        self.df_ref_trace = None

        self.region_qry = interval.region_qry

        self.pos = interval.region_ref.pos
        self.end = interval.region_ref.end

    def complete_anno_impl(self):

        # Get reference trace
        if self.df_ref_trace is None:
            self.df_ref_trace = get_reference_trace(self.interval)

        self.struct_qry = get_qry_struct_str(self.interval.df_segment)
        self.struct_ref = get_ref_struct_str(self.df_ref_trace)

    def row_impl(self):
        return pd.Series(
            [
                self.chrom, self.pos, self.end, self.variant_id, self.svtype, self.svlen,
                str(self.region_qry), self.strand,
                self.seg_n, self.struct_ref, self.struct_qry,
                self.score_variant,
                self.anchor_score_min, self.anchor_score_max,
                f'{self.chain_start_index}-{self.chain_end_index}'
            ],
            index=[
                '#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN',
                'QRY_REGION', 'QRY_STRAND',
                'SEG_N', 'STRUCT_REF', 'STRUCT_QRY',
                'VAR_SCORE',
                'ANCHOR_SCORE_MIN', 'ANCHOR_SCORE_MAX',
                'INTERVAL'
            ]
        )


class NullVariant(Variant):
    """
    Represents no variant call.
    """

    def __init__(self):
        Variant.__init__(self, None)

    def row_impl(self):
        raise RuntimeError('Null variant does not have a row')

class PatchVariant(Variant):
    """
    Represents no-variant call around alignment errors that do not represent real variation. For example, alignment
    artifacts around inverted repeats often break the alignment into pieces, but inspection of the pieces finds that
    they are largely contiguous with the reference. This patch variant type has a score of 0 and is treated as a bridge
    from one alignment anchor to another by graph traversal and prevents it from being pushed down non-optimal paths.
    """

    def __init__(self, interval=None):
        Variant.__init__(self, interval)

        self.score_variant = 0
        self.is_patch = True

        self.vartype = 'PATCH'


    def row_impl(self):
        raise RuntimeError('Null variant does not have a row')


class DevVariant(Variant):
    """
    A convenience class for development purposes. Not used by PAV.
    """

    def __init__(self, interval, caller_resources=None, var_region_kde=None):
        Variant.__init__(self, interval)

        self.vartype = 'DEV'
        self.caller_resources = caller_resources
        self.var_region_kde = var_region_kde

    def row_impl(self):
        raise RuntimeError('Dev variant does not have a row')


class VariantMatchError(Exception):
    """
    Raised by a Variant subclass when a variant call cannot be matched to a locus. Distinguishes "no variant call" from
    errors that might be raised by the constructor.
    """
    pass


#
# Reference and query structure routines
#

def get_reference_trace(interval):
    """
    Get a table representing the trace across the reference locus for this SV. This only covers the locus of the SV
    and omits distal template switches.

    In some cases, a complex SV has no reference context. If there is a clean insertion site (no deleted or duplciated
    bases at the SV site), but the SV is templated from multiple locations (or contains templated and untemplated
    insertions). In these cases, an empty table is returned.

    :param interval: Interval to generate a reference trace for.

    :return: A reference trace table. Has no rows if there is no reference context at the SV site.
    """

    # Set local template switch boundaries
    # Distance from reference positions (pos & end) must be no more half of:
    #   * CPX event in contig bases.
    #   * Distance between reference anchoring breakpoints.
    # May not expand beyond archoring alignments.
    local_dist = max([len(interval.region_ref), len(interval.region_qry)]) // 2

    local_pos = max([
        interval.region_ref.pos - local_dist,
        min([
            interval.df_segment.iloc[0]['POS'],
            interval.df_segment.iloc[0]['END'],
            interval.df_segment.iloc[-1]['POS'],
            interval.df_segment.iloc[-1]['END']
        ])
    ])

    local_end = min([
        interval.region_ref.end + local_dist,
        max([
            interval.df_segment.iloc[0]['POS'],
            interval.df_segment.iloc[0]['END'],
            interval.df_segment.iloc[-1]['POS'],
            interval.df_segment.iloc[-1]['END']
        ])
    ])

    # Make depth table
    df_depth = pavlib.align.util.align_bed_to_depth_bed(interval.df_segment.loc[interval.df_segment['IS_ALIGNED']], df_fai=None)

    # Refine depth to local regions (remove distal, trim segments crossing local region pos & end)
    df_depth['POS'] = df_depth['POS'].apply(lambda val: max(val, local_pos))
    df_depth['END'] = df_depth['END'].apply(lambda val: min(val, local_end))

    df_depth = df_depth.loc[
        (df_depth['#CHROM'] == interval.region_ref.chrom) &
        (df_depth['END'] - df_depth['POS'] > 0)
    ].copy()

    df_depth['INDEX'] = df_depth['INDEX'].fillna('').astype(str).apply(lambda val: val.strip())

    # Remove depth records that cover only anchor alignments at the flanks
    anchor_index_set = {
        str(interval.df_segment.iloc[0]['INDEX']),
        str(interval.df_segment.iloc[-1]['INDEX'])
    }

    while df_depth.shape[0] > 0 and df_depth.iloc[-1]['INDEX'] in anchor_index_set and df_depth.iloc[-1]['DEPTH'] == 1:
        df_depth = df_depth.iloc[:-1]

    while df_depth.shape[0] > 0 and df_depth.iloc[0]['INDEX'] in anchor_index_set and df_depth.iloc[0]['DEPTH'] == 1:
        df_depth = df_depth.iloc[1:]

    # QRY_ID and INDEX to tuples
    # df_depth['QRY_INDEX'] = df_depth.apply(lambda row:
    #     sorted(zip(
    #         row['QRY_ID'].split(','), row['INDEX'].split(',')
    #     )), axis=1
    # )
    #
    # del df_depth['QRY_ID']
    # del df_depth['INDEX']

    # Add forward and reverse counts
    df_depth['FWD_COUNT'] = 0
    df_depth['REV_COUNT'] = 0

    for index, row in df_depth.iterrows():

        depth = row['DEPTH']

        if depth > 0:
            index_set = {int(val) for val in row['INDEX'].split(',')}

            if len(index_set) != depth:
                raise RuntimeError(f'Depth record index list length {len(index_set)} does not match depth {depth}')

            fwd_count = np.sum(interval.df_segment.loc[interval.df_segment['INDEX'].isin(index_set), 'STRAND'] == '+')
            rev_count = depth - fwd_count

        else:
            fwd_count = 0
            rev_count = 0

        if interval.is_rev:
            fwd_count, rev_count = rev_count, fwd_count

        df_depth.loc[index, 'FWD_COUNT'] = fwd_count
        df_depth.loc[index, 'REV_COUNT'] = rev_count

    # # Aggregate records
    # if aggregate and df_depth.shape[0] > 0:
    #     df_depth_list = list()
    #     n = df_depth.shape[0]
    #
    #     i = 0
    #     row = df_depth.iloc[i]
    #
    #     while i < n:
    #
    #         if i + 2 < n:
    #             row2 = df_depth.iloc[i + 2]
    #
    #             if (
    #                 row['#CHROM'] == row2['#CHROM']
    #             ) and (
    #                 row['FWD_COUNT'] = row2
    #             ):
    #
    #     depth_len = df_depth['END'] - df_depth['POS']
    #
    #     for i in range(df_depth.shape[0] - 2):
    #
    #         if (
    #                 depth_len.iloc[i + 1] / (depth_len.iloc[i] + depth_len.iloc[i + 2])
    #         ) and (
    #             np.all(df_depth.iloc[i][['FWD_COUNT', 'REV_COUNT']] == df_depth.iloc[i + 2][['FWD_COUNT', 'REV_COUNT']])
    #         ):

    # Make reference context
    df_ref_trace_list = list()  # Type, length, depth fwd, depth rev, alignment index

    del_len = 0

    for index, row in df_depth.iterrows():

        depth_len = row['END'] - row['POS']
        depth = row['DEPTH']
        fwd_count = row['FWD_COUNT']
        rev_count = row['REV_COUNT']

        # Get segment type
        if depth == 0:
            seg_type = 'DEL'
            del_len += depth_len

        elif depth == 1:
            if rev_count == 1:
                seg_type = 'INV'
            else:
                seg_type = 'NML'

        else:
            dup_type = {
                0: 'DEL',
                1: 'NML',
                2: 'DUP',
                3: 'TRP',
                4: 'QUAD',
            }.get(depth, 'HDUP')  # Default to high-copy duplication

            if fwd_count == depth:
                seg_type = dup_type
            elif rev_count >= (depth - 1):
                seg_type = f'INV{dup_type}'
            else:
                seg_type = f'MIX{dup_type}'

        # Add trace record
        df_ref_trace_list.append(
            pd.concat([
                row,
                pd.Series(
                    [seg_type, depth_len],
                    index=['TYPE', 'LEN']
                )
            ], axis=0)
        )

    if len(df_ref_trace_list) == 0:
        return pd.DataFrame(
            [],
            columns=REF_TRACE_COLUMNS
        )

    df_ref_trace = pd.concat(df_ref_trace_list, axis=1).T

    if list(df_ref_trace.columns) != REF_TRACE_COLUMNS:
        raise RuntimeError(f'Unexpected reference trace columns (Program bug): {", ".join(df_ref_trace.columns)}')

    return df_ref_trace


def get_qry_struct_str(df_segment):
    """
    Get a comuplex SV structure following the reference through template switches, templated insertions, and
    untemplated insertions.

    :param df_segment: Segment table.

    :return: String describing the complex SV structure.
    """

    last_chrom = df_segment.iloc[0]['#CHROM']
    last_pos = df_segment.iloc[0]['END']

    struct_list = list()

    for i in range(1, df_segment.shape[0] - 1):
        row = df_segment.iloc[i]

        if row['IS_ALIGNED']:
            if row['#CHROM'] == last_chrom:
                dist = row['POS'] - last_pos
                struct_list.append(f'TS({dist})')
            else:
                struct_list.append(f'TS({row["#CHROM"]})')
                last_chrom = row['#CHROM']

            struct_list.append(f'TINS({row["STRAND"]}{row["LEN_QRY"]})')

            last_pos = row['END']
        else:
            struct_list.append(f'INS({row["LEN_QRY"]})')

    row = df_segment.iloc[-1]
    dist = row['POS'] - last_pos
    struct_list.append(f'TS({dist})')

    return ':'.join(struct_list)


def get_ref_struct_str(df_ref_trace):
    """
    Get reference structure string describing a complex SV from the reference perspective.

    :param df_ref_trace: Reference trace table.

    :return: A string describing the reference structure.
    """

    if df_ref_trace.shape[0] == 0:
        return 'INS'

    return ':'.join(df_ref_trace.apply(lambda row: f'{row["TYPE"]}({row["LEN"]})', axis=1))
