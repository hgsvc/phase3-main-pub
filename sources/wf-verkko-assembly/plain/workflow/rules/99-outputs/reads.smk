
READ_STATS = []

if RUN_COMPUTE_ALL_READ_STATS:
    READ_STATS.extend(
        rules.compute_all_read_stats.input.stats
    )
