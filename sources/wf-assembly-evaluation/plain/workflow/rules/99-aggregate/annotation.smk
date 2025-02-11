
ANNOTATION_OUTPUT = []

if RUN_HMMER:

    ANNOTATION_OUTPUT.extend(
        rules.run_all_hmmer_jobs.input.tables
    )

    ANNOTATION_OUTPUT.extend(
        rules.run_all_hmmer_jobs.input.stats
    )

    ANNOTATION_OUTPUT.extend(
        rules.run_all_hmmer_jobs.input.bed_out
    )


if True:
    # checking for N-gaps is a quasi zero-overhead
    # rule, i.e., not much of a point making this
    # a switch
    ANNOTATION_OUTPUT.extend(
        rules.run_all_ngaps_annotation.input.beds
    )

if RUN_REPEATMASKER:

    ANNOTATION_OUTPUT.extend(
        rules.run_all_repeatmasker_jobs.input.norm_tables
    )

    ANNOTATION_OUTPUT.extend(
        rules.run_all_repeatmasker_jobs.input.tars
    )


if True: # TODO: make proper switch RUN_SEQTK_TELOMERE:

    ANNOTATION_OUTPUT.extend(
        rules.run_all_seqtk_telomere_motifs.input.tables
    )
