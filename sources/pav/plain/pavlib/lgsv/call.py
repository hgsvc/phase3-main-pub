# Variant calling

import collections
import numpy as np
import os

import pavlib
import svpoplib


class VarRegionKde:
    def __init__(self, interval, caller_resources):

        # Default values
        self.df_kde = None
        self.df_rl = None
        self.kde_inv = False
        self.try_var = True

        if (
            interval.len_ref > 0
        ) and (
            np.min([interval.len_qry, interval.len_ref]) / np.max([interval.len_qry, interval.len_ref]) > 0.85
        ) and (
            interval.seg_n < 10
        ) and (
            interval.len_qry < 1e6
        ):

            # Match complex sequence to the gap region
            self.df_kde = pavlib.inv.get_state_table(
                region_ref=interval.region_ref,
                region_qry=interval.region_qry,
                ref_fa_name=caller_resources.ref_fa_name,
                qry_fa_name=caller_resources.qry_fa_name,
                ref_fai=caller_resources.ref_fai,
                qry_fai=caller_resources.qry_fai,
                is_rev=interval.region_qry.is_rev,
                k_util=caller_resources.k_util,
                kde=caller_resources.kde,
                max_ref_kmer_count=caller_resources.config_params.inv_max_ref_kmer_count,
                expand_bound=True,
                log=caller_resources.log_file
            )

            # Get run-length encoded table if at least half the k-mers were not discarded
            kmer_n = len(interval.region_qry) - caller_resources.k_util.k_size + 1

            if kmer_n > 0 and self.df_kde is not None and self.df_kde.shape[0] / kmer_n > 0.5:
                self.df_rl = pavlib.kde.rl_encoder(self.df_kde)

                if not np.any(self.df_rl['STATE'] == pavlib.inv.KDE_STATE_REV):
                    self.df_rl = None

            else:
                self.df_rl = None

            # Test for inverted significance
            if self.df_rl is not None:
                p_binom = pavlib.inv.test_kde(self.df_rl)  # Test KDE for inverted state significance

                # Test states for inverted or reference states
                qry_len = np.sum(self.df_rl['LEN_QRY'])

                prop_rev = np.sum(self.df_rl.loc[self.df_rl['STATE'] == pavlib.inv.KDE_STATE_REV, 'LEN_KDE']) / qry_len
                prop_fwdrev = np.sum(self.df_rl.loc[self.df_rl['STATE'] == pavlib.inv.KDE_STATE_FWDREV, 'LEN_KDE']) / qry_len
                prop_fwd = np.sum(self.df_rl.loc[self.df_rl['STATE'] == pavlib.inv.KDE_STATE_FWD, 'LEN_KDE']) / qry_len

                self.kde_inv = p_binom < 0.01 and prop_rev >= 0.5

                self.try_var = (
                    self.kde_inv
                ) or (
                    qry_len / kmer_n < 0.8  # Too many missing k-mers to make a call, try a variant
                ) or (
                    prop_fwd < 0.90 and (  # Too many forward states
                        prop_fwdrev < 0.5 or prop_fwd + prop_fwdrev < 0.90  # Weed out misalignments around inverted repeats (human chrX)
                    )
                )
            else:
                self.kde_inv = False

            # Check all segments, should belong to the gap region
            if self.kde_inv:
                qry_start = interval.df_segment.iloc[-1 if interval.is_rev else 0]['QRY_END']

                for index, row in interval.df_segment.loc[~ interval.df_segment['IS_ANCHOR']].iterrows():
                    pos, end = sorted([abs(row['QRY_POS'] - qry_start), abs(row['QRY_END'] - row['QRY_POS'])])

                    df_kde_seg = self.df_kde.loc[
                        (self.df_kde['INDEX'] >= pos) * (self.df_kde['INDEX'] <= end)
                    ]

                    df_kde_seg = df_kde_seg.loc[df_kde_seg['STATE'].isin({pavlib.inv.KDE_STATE_FWDREV, pavlib.inv.KDE_STATE_REV})]

                    prop_mer = df_kde_seg.shape[0] / (end - pos - caller_resources.k_util.k_size + 1)

                    if prop_mer < 0.5:
                        self.kde_inv = False


def call_from_align(caller_resources, min_anchor_score=pavlib.const.DEFAULT_MIN_ANCHOR_SCORE, verbose=False, dot_dirname=None):
    """
    Create a list of variant calls from alignment table.

    :param caller_resources: Caller resources.
    :param min_anchor_score: Minimum allowed score for an alignment segment to anchor a variant call.
    :param verbose: Print progress.
    :param dot_dirname: Directory where graph dot files are written.

    :return: A list of variant call objects.
    """

    variant_call_list = list()
    variant_id_set = set()

    min_anchor_score = pavlib.lgsv.util.get_min_anchor_score(min_anchor_score, caller_resources.score_model)

    for query_id in caller_resources.df_align_qry['QRY_ID'].unique():
        if verbose:
            print(f'Query: {query_id}')

        df_align = caller_resources.df_align_qry.loc[
            caller_resources.df_align_qry['QRY_ID'] == query_id
        ].reset_index(drop=True)

        if list(df_align['QRY_ORDER']) != list(df_align.index):
            raise RuntimeError(f'Query {query_id} order is out of sync with QRY_ORDER (Program Bug)')

        chain_set = set()

        start_index = 0
        last_index = df_align.shape[0]

        qryref_index_set = set(caller_resources.df_align_qryref.index)

        # Traverse interval setting each position to the left-most anchor candidate.
        while start_index < last_index:

            # Skip if anchor did not pass TIG & REF trimming
            if df_align.loc[start_index]['INDEX'] not in qryref_index_set:
                start_index += 1
                continue

            start_row = df_align.loc[start_index]
            end_index = start_index + 1

            # Traverse each interval after the start for right-most anchor candidates. Limit search by contig distance
            while end_index < last_index and pavlib.lgsv.chain.can_reach_anchor(start_row, df_align.loc[end_index], caller_resources.score_model):

                if df_align.loc[end_index]['INDEX'] in qryref_index_set and pavlib.lgsv.chain.can_anchor(
                        start_row, df_align.loc[end_index], caller_resources.score_model, min_anchor_score
                ):
                    chain_set.add((start_index, end_index))

                end_index += 1

            start_index += 1


        #
        # Resolve variant calls
        #

        sv_dict = dict()  # Key: interval range (tuple), value=SV object

        for start_index, end_index in chain_set:
            chain_node = pavlib.lgsv.chain.AnchorChainNode(start_index, end_index)

            interval = pavlib.lgsv.interval.AnchoredInterval(chain_node, df_align, caller_resources.score_model)

            # Get variant region KDE
            var_region_kde = VarRegionKde(interval, caller_resources)

            # Initialize to Null variant
            variant_call = pavlib.lgsv.variant.NullVariant()

            if var_region_kde.try_var:

                # Try INS
                variant_call = pavlib.lgsv.variant.try_variant(
                    pavlib.lgsv.variant.InsertionVariant, interval, caller_resources, variant_call, var_region_kde
                )

                # Try DEL
                variant_call = pavlib.lgsv.variant.try_variant(
                    pavlib.lgsv.variant.DeletionVariant, interval, caller_resources, variant_call, var_region_kde
                )

                # Try tandem duplication
                variant_call = pavlib.lgsv.variant.try_variant(
                    pavlib.lgsv.variant.TandemDuplicationVariant, interval, caller_resources, variant_call, var_region_kde
                )

                # Try Inversion
                variant_call = pavlib.lgsv.variant.try_variant(
                    pavlib.lgsv.variant.InversionVariant, interval, caller_resources, variant_call, var_region_kde
                )

                # Try complex
                variant_call = pavlib.lgsv.variant.try_variant(
                    pavlib.lgsv.variant.ComplexVariant, interval, caller_resources, variant_call, var_region_kde
                )
            else:
                variant_call = pavlib.lgsv.variant.PatchVariant(interval)

            # Set to Null variant if anchors cannot support the variant call
            if not variant_call.is_null() and variant_call.score_variant + variant_call.anchor_score_min < 0:
                variant_call = pavlib.lgsv.variant.NullVariant()

            # Complete candidate variant call
            sv_dict[(start_index, end_index)] = variant_call


        #
        # Choose optimal SVs
        #

        #df_align = caller_resources.df_align_qry

        # Initialize Bellman-Ford
        top_score = np.full(df_align.shape[0], -np.inf)   # Score, top-sorted graph
        top_tb = np.full(df_align.shape[0], -2)           # Traceback (points to parent node with the best score), top-sorted graph

        top_score[0] = 0
        top_tb[0] = -1

        start_index_list = list()

        # for source_node in chain_container.source_node_list:
        #     top_score[source_node.start_index] = 0
        #     top_tb[source_node.start_index] = -1  # Points to root node

        # if np.isneginf(top_score[0]):
        #     #raise RuntimeError('Root node does not point to node 0')
        #     print('Root node does not point to node 0')
        #     continue

        # Create a graph by nodes (anchor graph nodes are scored edges)
        node_link = collections.defaultdict(set)

        for start_index, end_index in chain_set:
            node_link[start_index].add(end_index)

        for start_index in range(df_align.shape[0] - 1):
            node_link[start_index].add(start_index + 1)

        # for start_index, end_index in sv_dict.keys():
        #     if not (sv_dict[start_index, end_index].is_null() or sv_dict[start_index, end_index].is_patch):
        #         node_link[start_index].add(end_index)

        # Update score by Bellman-Ford
        for start_index in range(df_align.shape[0]):
            base_score = top_score[start_index]

            # print(f'Start: {start_index:d}: {base_score:.2E}')  # DBGTMP

            if np.isneginf(base_score):  # Unreachable
                raise RuntimeError(f'Unreachable node at index {start_index}')

            for end_index in sorted(node_link[start_index]):

                # Score for this edge (Initialize to optimal score leading up to the edge)
                score = base_score

                sv_score = sv_dict[start_index, end_index].score_variant \
                    if (start_index, end_index) in sv_dict else -np.inf

                if not np.isneginf(sv_score):
                    # Variant call (or patch variant across alignment artifacts)
                    # Edge weight: Variant score and half of each anchor.
                    score +=  \
                        sv_score + \
                        df_align.loc[end_index]['SCORE'] / 2 + \
                        df_align.loc[start_index]['SCORE'] / 2

                else:
                    # No variant call
                    # Can be from in-chain edges (anchor candidates) or not-in-chain (edges added between sequential alignment records)
                    # Edge weight: Score by complex structure (template switches and gap penalties for each segment).

                    # Initial score: Template switches and gap penalties
                    score += \
                        pavlib.lgsv.variant.score_segment_transitions(
                            pavlib.lgsv.interval.get_segment_table(start_index, end_index, df_align, caller_resources.score_model),
                            caller_resources
                        )

                    if (start_index, end_index) in chain_set:
                        # In-chain edge
                        # Score: add half of anchor alignment scores and penalize by the reference gap
                        score += \
                            df_align.loc[end_index]['SCORE'] / 2 + \
                            df_align.loc[start_index]['SCORE'] / 2

                # print(f'\t* End: {end_index:d}: {score:.2E}')  # DBGTMP

                if score > top_score[end_index]:
                    # print('\t\t* New top')  # DBGTMP
                    top_score[end_index] = score
                    top_tb[end_index] = start_index

            # print() # DBGTMP
            # print('\n'.join([f'{r:2d}: {t:d} ({s:.2E})' for s, t, r in zip(top_score, top_tb, range(df_align.shape[0]))]))  # DBGTMP


        # # Trace back through scored paths
        # sink_node_list = sorted({
        #     chain_node.end_index
        #         for chain_node in chain_container.chain_dict.values()
        #             if len(chain_node.next_node_list) == 0
        # })
        #
        # if len(sink_node_list) == 0:
        #     continue

        last_node = df_align.shape[0] - 1

        # max_score_index = np.argmax([top_score[i] for i in sink_node_list])
        # max_sink_node = sink_node_list[max_score_index]
        #
        optimal_interval_list = list()
        #
        # last_node = max_sink_node

        while True:
            first_node = top_tb[last_node]

            if first_node < 0:
                break

            optimal_interval_list.append((first_node, last_node))

            last_node = first_node

        # Add variants and version IDs
        new_var_list = [
            sv_dict[node_interval] for node_interval in optimal_interval_list \
                if node_interval in sv_dict and not sv_dict[node_interval].is_null() and not sv_dict[node_interval].is_patch
        ]

        for variant_call in new_var_list:
            variant_call.variant_id = svpoplib.variant.version_id_name(variant_call.variant_id, variant_id_set)
            variant_id_set.add(variant_call.variant_id)

        variant_call_list.extend(new_var_list)

        # Write dot file
        if dot_dirname is not None:
            dot_filename = os.path.join(dot_dirname, f'lgsv_graph_{query_id}.dot.gz')

            with pavlib.io.PlainOrGzFile(dot_filename, 'wt') as out_file:
                pavlib.lgsv.util.dot_graph_writer(
                    out_file, df_align, chain_set, optimal_interval_list, sv_dict, graph_name=f'"{query_id}"', forcelabels=True
                )

    # Return list
    return variant_call_list
