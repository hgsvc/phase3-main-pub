
CONTIG_REF_ALIGN_OUTPUT = []

if True:  # TODO: no switch exists yet

    CONTIG_REF_ALIGN_OUTPUT.extend(
        rules.run_minimap_contig_to_ref_alignments.input.paf
    )
    CONTIG_REF_ALIGN_OUTPUT.extend(
        rules.run_minimap_contig_to_ref_alignments.input.bams
    )


if ASSEMBLY_UNITS_SEX_SPECIFIC:  # TODO: make proper switch

    CONTIG_REF_ALIGN_OUTPUT.extend(
        rules.build_assembly_karyotype_summary.output.tsv
    )
