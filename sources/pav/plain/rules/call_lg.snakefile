"""
Call alignment-truncating events (large SVs).
"""

import os
import pandas as pd
import shutil
import tarfile

import kanapy
import pavlib

global REF_FA
global get_config

global expand
global shell
global temp
global ASM_TABLE

# Call all large SVs
localrules: call_lg_all

rule call_lg_all:
    input:
        bed=lambda wildcards: pavlib.pipeline.expand_pattern(
            'results/{asm_name}/lgsv/svindel_ins_{hap}.bed.gz', ASM_TABLE, config
        )


# Call alignment-truncating SVs.
rule call_lg_discover:
    input:
        bed_qry='results/{asm_name}/align/trim-qry/align_qry_{hap}.bed.gz',
        bed_qryref='results/{asm_name}/align/trim-qryref/align_qry_{hap}.bed.gz',
        bed_none='results/{asm_name}/align/trim-none/align_qry_{hap}.bed.gz',
        fa_qry='data/query/{asm_name}/query_{hap}.fa.gz',
        fai_qry='data/query/{asm_name}/query_{hap}.fa.gz.fai',
        fa_ref='data/ref/ref.fa.gz',
        fai_ref='data/ref/ref.fa.gz.fai'
    output:
        bed_ins='results/{asm_name}/lgsv/svindel_ins_{hap}.bed.gz',
        bed_del='results/{asm_name}/lgsv/svindel_del_{hap}.bed.gz',
        bed_inv='results/{asm_name}/lgsv/sv_inv_{hap}.bed.gz',
        bed_cpx='results/{asm_name}/lgsv/sv_cpx_{hap}.bed.gz',
        bed_cpx_seg='results/{asm_name}/lgsv/segment_cpx_{hap}.bed.gz',
        bed_cpx_ref='results/{asm_name}/lgsv/reftrace_cpx_{hap}.bed.gz',
        dot_tar='results/{asm_name}/lgsv/lgsv_graph_{asm_name}_{hap}.tar'
    log:
        log='log/{asm_name}/lgsv/lgsv_call_{hap}.log'
    run:

        config_params = pavlib.config.get_config_dict(wildcards.asm_name, config, ASM_TABLE)

        # Set graph file output
        dot_dirname = f'temp/{wildcards.asm_name}/lgsv/graph_{wildcards.hap}'
        os.makedirs(dot_dirname, exist_ok=True)

        # Get score model
        score_model = pavlib.align.score.get_score_model(config_params.align_score_model)

        # Get minimum anchor score
        min_anchor_score = pavlib.lgsv.util.get_min_anchor_score(config_params.min_anchor_score, score_model)

        # Read alignments - Trim QRY
        df_align_qry = pd.read_csv(
            input.bed_qry,
            sep='\t',
            dtype={'#CHROM': str, 'QRY_ID': str}
        )

        df_align_qry.sort_values(['QRY_ID', 'QRY_ORDER'], inplace=True)
        df_align_qry.reset_index(inplace=True, drop=True)

        # Read alignments - Trim QRY/REF
        df_align_qryref = pd.read_csv(
            input.bed_qryref,
            sep='\t',
            index_col='INDEX',
            dtype={'#CHROM': str, 'QRY_ID': str}
        )

        # Read alignments - Trim NONE
        df_align_none = pd.read_csv(
            input.bed_none,
            sep='\t',
            index_col='INDEX',
            dtype={'#CHROM': str, 'QRY_ID': str}
        )

        # Get KDE for inversions
        kde = pavlib.kde.KdeTruncNorm(
            config_params.inv_kde_bandwidth, config_params.inv_kde_trunc_z, config_params.inv_kde_func
        )

        with open(log.log, 'w') as log_file:

            # Set caller resources
            caller_resources = pavlib.lgsv.util.CallerResources(
                df_align_qry=df_align_qry,
                df_align_qryref=df_align_qryref,
                df_align_none=df_align_none,
                qry_fa_name=input.fa_qry,
                ref_fa_name=input.fa_ref,
                hap=wildcards.hap,
                score_model=score_model,
                k_util=kanapy.util.kmer.KmerUtil(config_params.inv_k_size),
                kde=kde,
                log_file=log_file,
                config_params=config_params
            )

            # Call
            lgsv_list = pavlib.lgsv.call.call_from_align(caller_resources, min_anchor_score=min_anchor_score, dot_dirname=dot_dirname)

        # Remove duplicate IDs (most likely never actually removes variants, but creates table tracking issues if it does)
        lgsv_list_dedup = list()
        lgsv_id_set = set()

        for var in lgsv_list:
            if var.variant_id not in lgsv_id_set:
                lgsv_list_dedup.append(var)
                lgsv_id_set.add(var.variant_id)

        lgsv_list = lgsv_list_dedup

        # Create tables
        df_list = {
            'INS': list(),
            'DEL': list(),
            'INV': list(),
            'CPX': list()
        }

        for var in lgsv_list:
            row = var.row()

            if row['SVTYPE'] not in df_list.keys():
                raise RuntimeError(f'Unexpected SVTYPE: "{row["SVTYPE"]}"')

            df_list[row['SVTYPE']].append(row)

        if len(df_list['INS']) > 0:
            df_ins = pd.concat(df_list['INS'], axis=1).T
        else:
            df_ins = pd.DataFrame([], columns=pavlib.lgsv.variant.InsertionVariant(None, None).row().index)

        if len(df_list['DEL']) > 0:
            df_del = pd.concat(df_list['DEL'], axis=1).T
        else:
            df_del = pd.DataFrame([], columns=pavlib.lgsv.variant.DeletionVariant(None, None).row().index)

        if len(df_list['INV']) > 0:
            df_inv = pd.concat(df_list['INV'], axis=1).T
        else:
            df_inv = pd.DataFrame([], columns=pavlib.lgsv.variant.InversionVariant(None, None).row().index)

        if len(df_list['CPX']) > 0:
            df_cpx = pd.concat(df_list['CPX'], axis=1).T
        else:
            df_cpx = pd.DataFrame([], columns=pavlib.lgsv.variant.ComplexVariant(None, None).row().index)

        df_ins.sort_values(['#CHROM', 'POS', 'END', 'ID', 'QRY_REGION'], inplace=True)
        df_del.sort_values(['#CHROM', 'POS', 'END', 'ID', 'QRY_REGION'], inplace=True)
        df_inv.sort_values(['#CHROM', 'POS', 'END', 'ID', 'QRY_REGION'], inplace=True)
        df_cpx.sort_values(['#CHROM', 'POS', 'END', 'ID', 'QRY_REGION'], inplace=True)

        # Write variant tables
        df_ins.to_csv(output.bed_ins, sep='\t', index=False, compression='gzip')
        df_del.to_csv(output.bed_del, sep='\t', index=False, compression='gzip')
        df_inv.to_csv(output.bed_inv, sep='\t', index=False, compression='gzip')
        df_cpx.to_csv(output.bed_cpx, sep='\t', index=False, compression='gzip')

        # Write segment and reference trace tables
        df_segment_list = list()
        df_reftrace_list = list()

        for var in lgsv_list:
            if var.svtype != 'CPX':
                continue

            df_segment = var.interval.df_segment.copy()
            df_segment.insert(3, 'ID', var.variant_id)

            df_reftrace = var.df_ref_trace.copy()
            df_reftrace.insert(3, 'ID', var.variant_id)

            df_segment_list.append(df_segment)
            df_reftrace_list.append(df_reftrace)

        if len(df_segment_list) > 0:
            df_segment = pd.concat(df_segment_list, axis=0)
        else:
            df_segment = pd.DataFrame([], columns=(
                pavlib.lgsv.interval.SEGMENT_TABLE_FIELDS[:3] + ['ID'] + pavlib.lgsv.interval.SEGMENT_TABLE_FIELDS[3:]
            ))

        if len(df_reftrace_list) > 0:
            df_reftrace = pd.concat(df_reftrace_list, axis=0)
        else:
            df_reftrace = pd.DataFrame([], columns=(
                pavlib.lgsv.variant.REF_TRACE_COLUMNS[:3] + ['ID'] + pavlib.lgsv.variant.REF_TRACE_COLUMNS[3:]
            ))

        df_segment.to_csv(output.bed_cpx_seg, sep='\t', index=False, compression='gzip')
        df_reftrace.to_csv(output.bed_cpx_ref, sep='\t', index=False, compression='gzip')

        # Compress graph dot files
        if dot_dirname is not None:
            with tarfile.open(output.dot_tar,'w') as tar_file:
                for file in os.listdir(dot_dirname):
                    tar_file.add(os.path.join(dot_dirname, file))

            shutil.rmtree(dot_dirname)
