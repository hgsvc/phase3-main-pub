
rule compute_merged_label_count_statistics:
    input:
        regions = rules.add_ngap_sizes.output.table
    output:
        dump = DIR_RES.joinpath(
            "asm_label_qc", "merge_tables", "by-sample",
            "{sample}", "{sample}.merged-issues.{span}.augmented.count-stats.pck"
        )
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    run:
        import pandas as pd
        import collections as col
        import pickle as pck
        import itertools as itt

        df = pd.read_csv(input.regions, sep="\t", header=0)

        combo_counter = col.Counter()
        combo_counter["total_regions_num"] = df.shape[0]
        combo_counter["total_regions_bp"] = (df["end"] - df["start"]).sum()

        assembly_length = df.drop_duplicates("seq", inplace=False)["seq_length"].sum()
        ngap_length = df.drop_duplicates("seq", inplace=False)["ngap_length"].sum()
        combo_counter["total_assembly_length"] = assembly_length
        combo_counter["total_ngap_length"] = ngap_length
        combo_counter["total_gapless_length"] = assembly_length - ngap_length

        def count_labels(row, counter):

            if row.labels == "no-labels":
                # skip - no annotation
                return
            region_length = row.end - row.start
            total_labels = row.labels.split(",")
            distinct_labels = sorted(set(l.split("::")[0] for l in total_labels))
            num_distinct = len(distinct_labels)

            recorded_labels = set()
            # in the following, note that total_labels can contain
            # the same label several times if it was merged several
            # times into a larger region; hence, processing state
            # for each label is tracked via recorded_labels
            for label_bp in total_labels:
                label, bp = label_bp.split("::")
                bases = int(bp)

                # total stats
                if label not in recorded_labels:
                    counter[(label, "unpaired", "total", "region_num")] += 1
                    counter[(label, "unpaired", "total", "region_bp")] += region_length
                counter[(label, "unpaired", "total", "label_bp")] += bases

                if num_distinct == 1:
                    # singleton stats
                    counter[(label, "unpaired", "single", "region_num")] += 1
                    counter[(label, "unpaired", "single", "region_bp")] += region_length
                    counter[(label, "unpaired", "single", "label_bp")] += bases
                else:
                    # merged stats
                    if label not in recorded_labels:
                        counter[(label, "unpaired", "merged", "region_num")] += 1
                        counter[(label, "unpaired", "merged", "region_bp")] += region_length
                    counter[(label, "unpaired", "merged", "label_bp")] += bases

                recorded_labels.add(label)

            if num_distinct > 1:
                # combination stats - merged equals total
                for (a,b) in itt.combinations(distinct_labels, 2):
                    counter[("pair", (a, b), "region_num")] += 1
                    counter[("pair", (a, b), "region_bp")] += region_length
                counter[("combination", tuple(distinct_labels), "region_num")] += 1
                counter[("combination", tuple(distinct_labels), "region_bp")] += region_length

            return

        _ = df.apply(count_labels, axis=1, args=(combo_counter,))
        with open(output.dump, "wb") as cache:
            _ = pck.dump(combo_counter, cache)
    # END OF RUN BLOCK


rule compute_association_label_annotation:
    input:
        regions = rules.add_ngap_sizes.output.table
    output:
        table = DIR_RES.joinpath(
            "asm_label_qc", "merge_tables", "by-sample",
            "{sample}", "{sample}.merged-issues.{span}.augmented.assoc-tests.tsv"
        )
    conda:
        DIR_ENVS.joinpath("assessem.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    params:
        script=DIR_SCRIPTS.joinpath("asm_label_qc", "label_assoc.py")
    shell:
        "{params.script} --augmented-regions {input.regions} --output {output.table}"


localrules: merge_label_annotation_associations
rule merge_label_annotation_associations:
    input:
        tables = expand(
            rules.compute_association_label_annotation.output.table,
            sample=SAMPLES,
            allow_missing=True
        )
    output:
        table = DIR_RES.joinpath(
            "asm_label_qc", "merge_tables",
            f"SAMPLES.{ASSEMBLER}.merged-issues.{{span}}.augmented.assoc-tests.tsv"
        )
    run:
        import pandas as pd
        import pathlib as pl
        concat = []
        for table in sorted(input.tables):
            df = pd.read_csv(table, sep="\t", header=0, dtype={
                "fet_pvalue": str,
                "adj_pvalue": str,
            })
            sample = pl.Path(table).name.split(".merged-issues.")[0]
            df.insert(0, "sample", sample)
            concat.append(df)
        concat = pd.concat(concat, axis=0, ignore_index=False)
        concat.to_csv(output.table, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


localrules: extract_error_flag_global_summary
rule extract_error_flag_global_summary:
    input:
        dumps = expand(
            rules.compute_merged_label_count_statistics.output.dump,
            sample=SAMPLES,
            allow_missing=True
        )
    output:
        table = DIR_RES.joinpath(
            "asm_label_qc", "merge_tables",
            f"SAMPLES.{ASSEMBLER}.merged-issues.{{span}}.global-uniq.tsv"
        )
    run:
        import pandas as pd
        import collections as col
        import pickle as pck

        def extract_per_label_stats(counts):
            data_rows = []
            for key, count in counts.items():
                if isinstance(key, str):
                    if "total" in key:
                        data_rows.append(
                            ("assembly", key, count)
                        )
                if len(key) != 4:
                    continue
                assert key[1] == "unpaired"
                merged_key = key[2] + "_" + key[3].replace("region", "regions")
                data_rows.append(
                    (key[0], merged_key, count)
                )
            df = pd.DataFrame.from_records(
                data_rows, columns=["label", "statistic", "count"]
            )
            return df

        concat = []
        for cache_file in sorted(input.dumps):
            sample = pl.Path(cache_file).name.split(".merged-issues.")[0]
            au = pl.Path(cache_file).name.split(".")[3]
            with open(cache_file, "rb") as dump:
                counts = pck.load(dump)
            df = extract_per_label_stats(counts)
            df.insert(0, "sample", sample)
            df.insert(1, "assembly", au)
            concat.append(df)
        concat = pd.concat(concat, axis=0, ignore_index=False)
        concat.sort_values(["sample", "label"], inplace=True)
        concat.reset_index(drop=True, inplace=True)

        reduced_rows = col.defaultdict(dict)
        flag_labels = set()
        samples = set()
        for (sample, label), counts in concat.groupby(["sample", "label"]):
            samples.add(sample)
            if label == "assembly":
                prefix = "ASSM"
                stat_columns  = ["total_regions_num", "total_regions_bp", "total_gapless_length"]
            else:
                prefix = label
                stat_columns = [
                    "single_label_bp", "merged_label_bp", "single_regions_num",
                    "total_regions_num", "total_label_bp"
                ]
                flag_labels.add(label)
            for stat in stat_columns:
                try:
                    value = counts.loc[counts["statistic"] == stat, "count"].values[0]
                except IndexError:
                    value = 0
                reduced_rows[sample][f"{prefix}_{stat}"] = value
        for sample in samples:
            total_regions_num = reduced_rows[sample]["ASSM_total_regions_num"]
            total_regions_bp = reduced_rows[sample]["ASSM_total_regions_bp"]
            assm_size = reduced_rows[sample]["ASSM_total_gapless_length"]
            total_flagged_pct = round(total_regions_bp / assm_size * 100, 3)
            reduced_rows[sample]["ASSM_total_merged_flagged_pct"] = total_flagged_pct
            if "no-ont" not in wildcards.span:
                ont_single_bp = reduced_rows[sample]["ISPCON_single_label_bp"]
                adj_flagged_bp = total_regions_bp - ont_single_bp
                adj_flagged_pct = round(adj_flagged_bp / assm_size * 100, 3)
                reduced_rows[sample]["ASSM_total_merged_flagged_ontfree_pct"] = adj_flagged_pct

            for label in flag_labels:
                try:
                    label_single_regions = reduced_rows[sample][f"{label}_single_regions_num"]
                except KeyError:  # happens for SSQBRK
                    uniq_pct = 0.
                else:
                    uniq_pct = round(label_single_regions / total_regions_num * 100, 3)
                reduced_rows[sample][f"{label}_regions_uniq_pct"] = uniq_pct

        df = pd.DataFrame.from_dict(reduced_rows, orient="index")
        df.fillna(0., inplace=True)

        df.to_csv(output.table, sep="\t", header=True, index=True, index_label="sample")
    # END OF RUN BLOCK


rule run_all_merged_region_stats:
    input:
        dumps = expand(
            rules.compute_merged_label_count_statistics.output.dump,
            sample=SAMPLES,
            span=["wg", "ps", "wg-no-ont", "ps-no-ont"]
        ),
        tables = expand(
            rules.compute_association_label_annotation.output.table,
            sample=SAMPLES,
            span=["wg", "ps", "wg-no-ont", "ps-no-ont"]
        ),
        merged_assoc = expand(
            rules.merge_label_annotation_associations.output.table,
            span=["wg", "ps", "wg-no-ont", "ps-no-ont"]
        ),
        merged_global = expand(
            rules.extract_error_flag_global_summary.output.table,
            span=["wg", "ps", "wg-no-ont", "ps-no-ont"]
        )
