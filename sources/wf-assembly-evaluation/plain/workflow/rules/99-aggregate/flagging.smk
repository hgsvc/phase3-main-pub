
FLAG_REGIONS_OUTPUT = []

# run nucfreq
if True:  # TODO: make proper switch

    FLAG_REGIONS_OUTPUT.extend(
        rules.run_all_nucfreq_jobs.input.hdf
    )
    FLAG_REGIONS_OUTPUT.extend(
        rules.run_all_nucfreq_jobs.input.regions
    )

    FLAG_REGIONS_OUTPUT.extend(
        rules.run_all_deepvariant_hifi_mismatches.input.vcf
    )
