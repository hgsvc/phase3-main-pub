
rule prepare_window_read_coverage_histogram:
    input:
        hdf = rules.transform_mosdepth_window_read_coverage.output.hdf
    output:
        tsv = DIR_PROC.joinpath(
            "60-flagging", "readcov", "window_histogram",
            "{sample}.read-cov-hist.tsv"
        )
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    params:
        low_bin=-5,
        high_bin=305,
        step_size=5,
        mapq=0,
        signal="pct_median_cov"
    run:
        import numpy as np
        import pandas as pd
        import pathlib as pl
        import sys

        bins = np.arange(params.low_bin, params.high_bin, step=params.step_size).round(0)
        bins = np.hstack((bins, np.array([params.high_bin, sys.maxsize])))

        bin_labels = []
        for edge_idx, bin_edge in enumerate(bins[1:], start=1):
            if edge_idx == bins.size - 1:
                label = "LRG"
            else:
                label = f"{bin_edge:03}"
            bin_labels.append(label)

        # TODO move to settings module
        cov_groups = GROUPS_WINDOW_READ_COVERAGE

        group_counts = dict()
        with pd.HDFStore(input.hdf, "r") as hdf:
            for key in hdf.keys():
                use_group = None
                for group_name in cov_groups:
                    if key.startswith(f"/{group_name}"):
                        use_group = group_name
                        break
                if use_group is None:
                    continue

                data = hdf[key]
                for column in data.columns:
                    read_type, aln_type, mapq, stat = column
                    if stat != params.signal:
                        continue
                    if mapq != params.mapq:
                        continue

                    hist = pd.cut(
                        data.loc[:, column],
                        bins=bins,
                        include_lowest=False, right=True,
                        labels=bin_labels
                    ).value_counts()
                    store_key = (wildcards.sample, read_type, aln_type, mapq, use_group)
                    if store_key not in group_counts:
                        tmp = pd.Series(
                            np.zeros(bins.size-1, dtype=int),
                            index=bin_labels
                        )
                        group_counts[store_key] = tmp
                    group_counts[store_key] += hist

        sample_hist = pd.DataFrame.from_dict(group_counts, orient="index")
        sample_hist.index.rename(["sample", "read_type", "aln_type", "mapq", "asm_unit"], inplace=True)
        sample_hist.to_csv(output.tsv, sep="\t", header=True, index=True)
    # END OF RUN BLOCK


rule extract_issue_windows_readcov_empty:
    input:
        hdf = rules.transform_mosdepth_window_read_coverage.output.hdf
    output:
        bed = DIR_RES.joinpath(
            "regions", "{sample}",
            "{sample}.no-read-support.bed.gz"
        )
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt
    params:
        low_cov_all = 1
    run:
        import pandas as pd
        import numpy as np
        import gzip as gz

        cov_groups = GROUPS_WINDOW_READ_COVERAGE
        read_types = USE_READ_TYPES_WINDOW_COVERAGE

        dump_regions = []

        with pd.HDFStore(input.hdf, "r") as hdf:
            for key in hdf.keys():
                has_group = None
                for cov_group in cov_groups:
                    if cov_group in key:
                        has_group = cov_group
                        break
                if has_group is None:
                    continue
                cov_data = hdf[key]

                empty_regions = cov_data.loc[(cov_data < params.low_cov_all).all(axis=1), :].copy()
                if not empty_regions.empty:
                    empty_regions["name"] = "no_read_support"
                    empty_regions.set_index([empty_regions.index, "name"], inplace=True)
                    dump_regions.append(empty_regions)

        if len(dump_regions) > 0:
            dump_regions = pd.concat(dump_regions, axis=0, ignore_index=False)
            dump_regions.sort_index(inplace=True)
            dump_regions.columns = flatten_window_readcov_columns(dump_regions.columns.names, dump_regions.columns)
            dump_regions.reset_index(drop=False, inplace=True)

            # by construction, the contigs are tagged;
            # turn back into regular contig names here
            dump_regions.rename({"contig": "tagged_contig"}, axis=1, inplace=True)
            # everything except 0/contig - 1/start - 2/end - 3/name
            value_columns = list(dump_regions.columns[4:])
            dump_regions["contig"] = dump_regions["tagged_contig"].apply(lambda x: x.rsplit(".", 1)[0])
            dump_regions["asm_unit"] = dump_regions["tagged_contig"].apply(lambda x: x.rsplit(".", 1)[1])
            dump_regions.drop("tagged_contig", axis=1, inplace=True)
            dump_regions = dump_regions[["contig", "start", "end", "name", "asm_unit"] + value_columns]

            with gz.open(output.bed, "wt") as bed:
                _ = bed.write("#")
                dump_regions.to_csv(bed, sep="\t", header=True, index=False)
        else:
            with gz.open(output.bed, "wt") as bed:
                _ = bed.write("#contig\tstart\tend\tname\n")

    # END OF RUN BLOCK


rule extract_issue_windows_readcov_onetype:
    input:
        hdf = rules.transform_mosdepth_window_read_coverage.output.hdf
    output:
        bed = DIR_RES.joinpath(
            "regions", "{sample}",
            "{sample}.only-{read_type}-support.bed.gz"
        )
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt
    params:
        min_cov_pct_reads = 10,
        max_cov_pct_other = 1
    run:
        import pandas as pd
        import numpy as np
        import gzip as gz

        cov_groups = GROUPS_WINDOW_READ_COVERAGE
        read_types = USE_READ_TYPES_WINDOW_COVERAGE
        assert wildcards.read_type in read_types

        dump_regions = []

        with pd.HDFStore(input.hdf, "r") as hdf:
            for key in hdf.keys():
                has_group = None
                for cov_group in cov_groups:
                    if cov_group in key:
                        has_group = cov_group
                        break
                if has_group is None:
                    continue
                cov_data = hdf[key]

                region_indexer = np.zeros(cov_data.shape[0], dtype=bool)
                region_indexer[:] = False
                supported_regions = ((cov_data.xs(
                    (wildcards.read_type, "pct_median_cov"),
                    level=["read_type", "value"],
                    axis=1
                )) > params.min_cov_pct_reads).any(axis=1).values
                # indexer: region has READ coverage
                region_indexer[supported_regions] = True

                for other_read_type in read_types:
                    if other_read_type == wildcards.read_type:
                        continue
                    dropout_regions = ((cov_data.xs(
                        (other_read_type, "pct_median_cov"),
                        level=["read_type", "value"],
                        axis=1
                    )) < params.max_cov_pct_other).all(axis=1).values
                    # indexer: region has READ coverage AND no OTHER coverage
                    region_indexer = np.logical_and(
                        region_indexer,
                        dropout_regions
                    )

                single_regions = cov_data.loc[region_indexer, :].copy()
                if not single_regions.empty:
                    single_regions["name"] = f"only_{wildcards.read_type}_support"
                    single_regions.set_index([single_regions.index, "name"], inplace=True)
                    dump_regions.append(single_regions)

        if len(dump_regions) > 0:
            dump_regions = pd.concat(dump_regions, axis=0, ignore_index=False)
            dump_regions.sort_index(inplace=True)
            dump_regions.columns = flatten_window_readcov_columns(dump_regions.columns.names, dump_regions.columns)
            dump_regions.reset_index(drop=False, inplace=True)

            # by construction, the contigs are tagged;
            # turn back into regular contig names here
            dump_regions.rename({"contig": "tagged_contig"}, axis=1, inplace=True)
            # everything except 0/contig - 1/start - 2/end - 3/name
            value_columns = list(dump_regions.columns[4:])
            dump_regions["contig"] = dump_regions["tagged_contig"].apply(lambda x: x.rsplit(".", 1)[0])
            dump_regions["asm_unit"] = dump_regions["tagged_contig"].apply(lambda x: x.rsplit(".", 1)[1])
            dump_regions.drop("tagged_contig", axis=1, inplace=True)
            dump_regions = dump_regions[["contig", "start", "end", "name", "asm_unit"] + value_columns]

            with gz.open(output.bed, "wt") as bed:
                _ = bed.write("#")
                dump_regions.to_csv(bed, sep="\t", header=True, index=False)
        else:
            with gz.open(output.bed, "wt") as bed:
                _ = bed.write("#contig\tstart\tend\tname\n")

    # END OF RUN BLOCK


rule run_all_window_read_coverage_histograms:
    input:
        tsv = expand(
            rules.prepare_window_read_coverage_histogram.output.tsv,
            sample=SAMPLES
        )


rule run_all_extract_window_readcov_issue_regions:
    input:
        empty = expand(
            rules.extract_issue_windows_readcov_empty.output.bed,
            sample=SAMPLES
        ),
        single = expand(
            rules.extract_issue_windows_readcov_onetype.output.bed,
            read_type=USE_READ_TYPES_WINDOW_COVERAGE,
            sample=SAMPLES
        )
