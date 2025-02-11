
ASSEMBLY_COMPLETENESS_OUTPUT = []

if True:  # TODO: make proper switch

    ASSEMBLY_COMPLETENESS_OUTPUT.extend(
        rules.run_all_asm_completeness_asmgene.input
    )

    ASSEMBLY_COMPLETENESS_OUTPUT.extend(
        rules.run_all_compute_approx_ref_span.input.tables
    )

    ASSEMBLY_COMPLETENESS_OUTPUT.extend(
        rules.run_all_compleasm.input
    )

    ASSEMBLY_COMPLETENESS_OUTPUT.extend(
        rules.run_all_label_contig_alignments.input
    )
