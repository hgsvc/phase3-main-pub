
CONTAMINATION_OUTPUT = []

if RUN_NCBI_FCS_ADAPTOR or RUN_NCBI_FCS_GX:
    # TODO
    # due to the filtering script, both outputs
    # have to exist. Running just the adaptor screening
    # seems anyway an unlikely thing to do, so the
    # implicit coupling here may not be harmful

    CONTAMINATION_OUTPUT.extend(
        rules.run_all_ncbi_fcs_reports.input.rep_adapter
    )
    CONTAMINATION_OUTPUT.extend(
        rules.run_all_ncbi_fcs_reports.input.rep_gxcontam
    )

    CONTAMINATION_OUTPUT.extend(
        rules.run_all_clean_assembly.input.asm_units
    )
    CONTAMINATION_OUTPUT.extend(
        rules.run_all_clean_assembly.input.contam
    )
    CONTAMINATION_OUTPUT.extend(
        rules.run_all_clean_assembly.input.regions_tagged
    )
    CONTAMINATION_OUTPUT.extend(
        rules.run_all_clean_assembly.input.regions_untagged
    )
