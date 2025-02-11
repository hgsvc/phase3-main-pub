
READ_ASSEMBLY_ALIGN_OUTPUT = []

if True:  # TODO: make switch - ? Only thing that works w/o reference

    READ_ASSEMBLY_ALIGN_OUTPUT.extend(
        rules.run_all_get_primary_alignment_read_lists.input.tsv
    )

    READ_ASSEMBLY_ALIGN_OUTPUT.extend(
        rules.run_all_window_read_coverage_histograms.input.tsv
    )

    READ_ASSEMBLY_ALIGN_OUTPUT.extend(
        rules.run_all_extract_window_readcov_issue_regions.input.empty
    )

    READ_ASSEMBLY_ALIGN_OUTPUT.extend(
        rules.run_all_extract_window_readcov_issue_regions.input.single
    )
