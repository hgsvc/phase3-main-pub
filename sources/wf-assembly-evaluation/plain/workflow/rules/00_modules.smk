"""
Use this module to list all includes
required for your pipeline - do not
add your pipeline-specific modules
to "commons/00_commons.smk"
"""

include: "00-prepare/settings.smk"
include: "00-prepare/sample_table.smk"
include: "00-prepare/source_setup.smk"

include: "05-preprocess/filter_unit.smk"
include: "05-preprocess/merge_tag.smk"
include: "05-preprocess/ncbi_fcs.smk"
include: "05-preprocess/remove_contam.smk"
include: "05-preprocess/ngaps.smk"

include: "10-asm-align/pyutils.smk"
include: "10-asm-align/align.smk"

include: "20-read-align/align.smk"
include: "20-read-align/filter.smk"
include: "20-read-align/merge.smk"

include: "50-postprocess/pyutils.smk"
include: "50-postprocess/asm_chrom_assign.smk"
include: "50-postprocess/asm_karyo_est.smk"
include: "50-postprocess/asm_ctg_refcov.smk"
include: "50-postprocess/asm_ctg_readcov.smk"

include: "60-flagging/pyutils.smk"
include: "60-flagging/nucfreq.smk"
include: "60-flagging/readcov.smk"
include: "60-flagging/misjoins.smk"
include: "60-flagging/mismatches.smk"

include: "70-annotate/pyutils.smk"
include: "70-annotate/setup_repmask.smk"
include: "70-annotate/repeatmasker.smk"
include: "70-annotate/hmmer.smk"
include: "70-annotate/telomeres.smk"

include: "75-completeness/pyutils.smk"
include: "75-completeness/asmgene.smk"
include: "75-completeness/busco.smk"
include: "75-completeness/ref_span.smk"
include: "75-completeness/breaks.smk"

include: "80-statistics/assemblies.smk"

include: "99-aggregate/annotation.smk"
include: "99-aggregate/completeness.smk"
include: "99-aggregate/contamination.smk"
include: "99-aggregate/contig_align.smk"
include: "99-aggregate/flagging.smk"
include: "99-aggregate/mosdepth.smk"
include: "99-aggregate/read_align.smk"
include: "99-aggregate/statistics.smk"
