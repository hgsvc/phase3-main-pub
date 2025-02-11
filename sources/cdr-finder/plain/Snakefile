import pandas as pd
import numpy as np

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

configfile: "config.yaml"

MANIFEST = config.get("manifest", "manifest.tab")
WINDOW = config["window_size"]
SLIDE = config["slide"]
REPORT_THRESHOLD = config.get("report_threshold", 25000)
ALR_THRESHOLD = config.get("alr_threshold", 100000)
EDGE_SEARCH = config.get("edge_search", 50000)
LOW_THRESHOLD = config.get("low_threshold", 39)


manifest_df = pd.read_csv(MANIFEST, sep='\t', index_col="SAMPLE", dtype=str)

def find_fasta(wildcards):
    return manifest_df.at[wildcards.sample, 'FASTA']

def find_target_bed(wildcards):
    return manifest_df.at[wildcards.sample, 'TARGET_BED']

def find_meth_tsv(wildcards):
    return manifest_df.at[wildcards.sample, 'METH_TSV']

def getOverlap(alin_rec):
    a = [alin_rec["start"], alin_rec["end"]]
    b = [alin_rec["bin_start"], alin_rec["bin_stop"]]
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

wildcard_constraints:
    sample="|".join(manifest_df.index)

rule all:
    input:
        expand("results/{sample}_CDR.bed", sample=manifest_df.index),


###
# Subset fasta: Extracts subsequence from target bed
###

rule subset_fasta:
    input:
        fasta=find_fasta,
        bed=find_target_bed,
    output:
        rm_fasta="results/fasta/{sample}_subset.fasta",
    resources:
        mem=8,
        hrs=1,
    threads: 1
    singularity:
        "docker://eichlerlab/binf-basics:0.1"
    shell:
        """
        bedtools getfasta -fi {input.fasta} -bed {input.bed} > {output.rm_fasta}
        """

###
# Subset meth: Intersects target bed with methylation bed
###


rule subset_meth:
    input:
        methylation_tsv=find_meth_tsv,
        target_bed=find_target_bed,
    output:
        subset_bed="results/{sample}_subset.bed",
    resources:
        mem=8,
        hrs=12,
    threads: 1
    singularity:
        "docker://eichlerlab/binf-basics:0.1"
    shell:
        """
        bedtools intersect -a {input.methylation_tsv} -b {input.target_bed} -wa > {output.subset_bed}
        """

###
# run_rm: Runs repeatmasker on the extracted subsequence fasta from rule subset_fasta
###


rule run_rm:
    input:
        fasta=rules.subset_fasta.output.rm_fasta,
    output:
        rm_out="results/rm/{sample}_subset.fasta.out",
    resources:
        mem=8,
        hrs=12,
    threads: 12
    singularity:
        "docker://eichlerlab/assembly_eval:0.1"
    shell:
        """
        RepeatMasker -species human -dir $( dirname {output.rm_out} ) -qq -pa {threads} {input.fasta}
        """

###
# Calculates mean frequency in windows/bins of the methylation tsv over the target region
###

rule calc_windows:
    input:
        methylation_tsv=rules.subset_meth.output.subset_bed,
        target_bed=find_target_bed,
    output:
        binned_freq="results/{sample}_binned_freq.bed",
    resources:
        mem=8,
        hrs=1,
    threads: 1
    params:
        window_size=WINDOW,
        slide=SLIDE,
    run:
        bed_df = pd.read_csv(
            input.target_bed, header=None, sep="\t", names=["chrom", "start", "stop"]
        )
        meth_df = pd.read_csv(
            input.methylation_tsv,
            header=None,
            sep="\t",
            names=[
                "chrom",
                "start",
                "end",
                "mod_base_name",
                "score",
                "strand",
                "n1",
                "n2",
                "n3",
                "coverage",
                "freq",
                "ncanon",
                "nmod",
                "nfilt",
                "nnocall",
                "naltmod",
            ],
        )
        keep_window = pd.DataFrame()

        for index in bed_df.index:
            start = bed_df.at[index, "start"]
            stop = bed_df.at[index, "stop"]
            chrom = bed_df.at[index, "chrom"]
            chrom_meth_df = meth_df.loc[meth_df["chrom"] == chrom].copy()
            seq = np.arange(
                int(start), int(stop) + int(params.slide), int(params.slide)
            )
            window_size = (int(int(params.window_size) / int(params.slide))) + 1
            window_avg = pd.DataFrame(columns=["Chr", "Bin", "Freq"])

            for i in range(len(seq) - window_size + 1):
                window = seq[i : i + window_size]
                window_meth_pre = chrom_meth_df.copy()
                RM_start = window[0] - start
                RM_stop = window[-1] - start
                window_meth_pre["bin_start"] = window[0]
                window_meth_pre["bin_stop"] = window[-1]
                print(RM_start, RM_stop)
                if len(window_meth_pre) == 0:
                    continue
                window_meth_pre["overlap"] = window_meth_pre.apply(getOverlap, axis=1)
                window_meth = window_meth_pre.loc[window_meth_pre["overlap"] != 0]
                if len(window_meth) == 0:
                    continue
                avg_meth = window_meth["freq"].mean()
                new_row = pd.DataFrame({
                        "Chr": [str(chrom)],
                        "Bin": [str(RM_start) + "-" + str(RM_stop)],
                        "Freq": [avg_meth]
                })
                window_avg = pd.concat([window_avg, new_row], ignore_index=True)

            window_avg = window_avg.round({"Freq": 2})
            keep_window = pd.concat([keep_window, window_avg])

        keep_window.reset_index(inplace=True)

        with open(output.binned_freq, "w+") as outfile:
            for idx_2 in keep_window.index:
                pos1 = str(keep_window.at[idx_2, "Bin"]).split("-")[0]
                pos2 = str(keep_window.at[idx_2, "Bin"]).split("-")[1]
                freq = str(keep_window.at[idx_2, "Freq"])
                chrom = str(keep_window.at[idx_2, "Chr"])
                outfile.write("\t".join([chrom, pos1, pos2, freq, f"{idx_2}\n"]))


###
# Format RM: Converts Repeatmasker output to bedfile and extracts only ALR annotations
###

rule format_RM:
    input:
        repeat_masker=rules.run_rm.output.rm_out,
    output:
        repeat_bed="results/{sample}_rm.bed",
    resources:
        mem=8,
        hrs=1,
    threads: 1
    shell:
        """
        awk -v OFS="\\t" '{{print $5, $6, $7, $10, $9}}' {input.repeat_masker} | grep "ALR" > {output.repeat_bed}
        """

###
# Filter RM: Returns Repeatmasker bed contig name to original fasta name
###


rule filter_RM:
    input:
        repeat_masker=rules.format_RM.output.repeat_bed,
    output:
        rm_bed="results/{sample}_rm_ALR.bed",
    resources:
        mem=8,
        hrs=1,
    threads: 1
    run:
        RM_df = pd.read_csv(
            input.repeat_masker,
            header=None,
            sep="\t",
            names=["chr", "start", "stop", "ALR", "or"],
        )
        RM_df["chr"] = RM_df["chr"].str.split(pat=":").str[0]
        RM_df.to_csv(output.rm_bed, header=None, index=None, sep="\t")

###
# Intersect RM: Merges repeatmasker bed file with 500bp slop, and only keeps those above the ALR Threshold
#               Intersects merged ALR bed file with binned methylation frequency bedfile 
###


rule intersect_RM:
    input:
        repeat_masker=rules.filter_RM.output.rm_bed,
        binned_freq=rules.calc_windows.output.binned_freq,
    output:
        intersect_bed="results/{sample}_intersect.bed",
        merged="results/{sample}_rm_merged.bed",
    resources:
        mem=8,
        hrs=1,
    params:
        alr_threshold=ALR_THRESHOLD
    threads: 1
    singularity:
        "docker://eichlerlab/binf-basics:0.1"
    shell:
        """
        bedtools merge -i {input.repeat_masker} -d 500 | awk '{{if ($3-$2 > {params.alr_threshold}) print}}' > {output.merged}
        bedtools intersect -a {input.binned_freq} -b {output.merged} -f 1 -wa -u > {output.intersect_bed}
        """

###
# in threshold: Evaluates bed file from intersected bed file and surrounding regions for signatures of CDR
# Initial signature is mean methylation frequency below dynamically calculated threshold 
# Confidence is evaluated based on whether or not these low methylation frequency windows clear a size threshold 
# + have flanking peaks of high methylation 
###


rule in_threshold:
    input:
        intersect_bed=rules.intersect_RM.output.intersect_bed,
        target_bed=find_target_bed,
    output:
        final_call="results/{sample}_CDR.bed",
        windows_call="results/{sample}_windows_CDR.bed",
    resources:
        mem=8,
        hrs=1,
    threads: 1
    run:
        # read in intersection bin bed
        intersect_bed = pd.read_csv(
            input.intersect_bed,
            header=None,
            sep="\t",
            names=["chrom", "start", "stop", "freq", "idx"],
            index_col="idx",
        )
        with open(input.target_bed) as infile, open(output.final_call, "w+") as outfile:
            # Evaluate over all windows in target bed        
            for line in infile:
                # Empty dataframe for windows
                window_avg_all = pd.DataFrame(columns=["Bin", "Freq"])
                # Parse target bed
                chrom, start = line.split("\t")[0], line.split("\t")[1]
                # Subset to chromosome
                intersect_bed_sub = intersect_bed.loc[intersect_bed["chrom"] == chrom]
                # Convert intersect bed to expected format for df
                for index in intersect_bed_sub.index:
                    chrom_int = intersect_bed_sub.at[index, "chrom"]
                    start_int = int(intersect_bed_sub.at[index, "start"]) + int(start)
                    stop_int = int(intersect_bed_sub.at[index, "stop"]) + int(start)
                    freq = intersect_bed_sub.at[index, "freq"]
                    new_row2 = pd.DataFrame({
                        "Bin": [str(start_int) + "-" + str(stop_int)],
                        "pos": [start_int],
                        "end": [stop_int],
                        "Freq": [float(freq)],
                        "idx_2": [index]
                    })
                    window_avg_all = pd.concat([window_avg_all, new_row2], ignore_index=True)

                window_avg_all["idx_2"] = window_avg_all["idx_2"].astype(int)
                window_avg_all = window_avg_all.set_index("idx_2")
                # Establish threshold, stddev, max
                median_freq = np.median(window_avg_all["Freq"])
                stddev_freq = np.std(window_avg_all["Freq"])
                max_freq = max(window_avg_all["Freq"])
                mean_freq = np.mean(window_avg_all["Freq"])
                max_thresh = max_freq - 1 * stddev_freq
                # Subset df to those windows under the threshold for instances where mean represents high methylation (large centromeres)
                if mean_freq > LOW_THRESHOLD:
                    threshold = mean_freq - 1 * stddev_freq
                    window_avg_under = window_avg_all.loc[
                        window_avg_all["Freq"] <= threshold
                    ].copy()
                # If methylation mean is low (small centromeres with more windows in CDR than not)
                else:
                    threshold = mean_freq + 1 * (stddev_freq/2)
                    window_avg_under = window_avg_all.loc[
                        ~(window_avg_all["Freq"] > threshold)
                    ].copy()
                s = pd.Series(window_avg_under.index.tolist())
                # Puts windows into consecutive ranges
                groups = (
                    s.groupby(s.diff().ne(1).cumsum())
                    .apply(
                        lambda x: [x.iloc[0], x.iloc[-1]]
                        if len(x) >= 2
                        else [x.iloc[0]]
                    )
                    .tolist()
                )
                groups = [sub for sub in groups if len(sub) > 1]
                # Begin checking for content of gaps
                check_df = pd.DataFrame()
                # Make df containing only windows which are consecutive and under threshold
                for y in groups:
                    check_df = pd.concat([check_df, window_avg_all.loc[y[0] : y[1]]])
                # get iterable of indexes
                initial_index_list = check_df.index
                # Begin checking for gaps in the ranges
                for i, index in enumerate(initial_index_list[:-1]):
                    # If next window found,continue
                    if index + 1 in initial_index_list:
                        continue
                    # If not, make df containing only the windows between the gaps
                    else:
                        gap_df = window_avg_all.loc[
                            initial_index_list[i] + 1 : initial_index_list[i + 1] - 1
                        ]
                        # If one window is within one stddev of the max, do not add gap
                        if (
                            len(gap_df.loc[gap_df["Freq"] > (max_freq - stddev_freq)])
                            == 0
                        ):
                            check_df = pd.concat([check_df, gap_df])
                        else:
                            continue
                s = pd.Series(sorted(check_df.index.tolist()))
                # Puts windows into consecutive ranges, again
                groups = (
                    s.groupby(s.diff().ne(1).cumsum())
                    .apply(
                        lambda x: [x.iloc[0], x.iloc[-1]]
                        if len(x) >= 2
                        else [x.iloc[0]]
                    )
                    .tolist()
                )
                groups = [sub for sub in groups if len(sub) > 1]
                # Writes output
                window_df_out = pd.DataFrame()
                for y in groups:
                    pos1 = int(check_df.at[y[0], "Bin"].split("-")[0])
                    pos2 = int(check_df.at[y[1], "Bin"].split("-")[1])
                    # left_check and right_check are checking for flanking peaks of high methylation frequency
                    left_check = window_avg_all.loc[(window_avg_all['end'] < pos1 ) & (window_avg_all['end'] >= pos1-EDGE_SEARCH )].copy()
                    right_check = window_avg_all.loc[(window_avg_all['end'] > pos2 ) & (window_avg_all['end'] <= pos2+EDGE_SEARCH )].copy()
                    # If neighbors are found on either side label as high confidence
                    if np.max(right_check['Freq']) > max_thresh and np.max(left_check['Freq']) > max_thresh:
                        confidence = "high_confidence"
                    # Else assign and label low confidence
                    else:
                        confidence = "low_confidence-neighbor_peaks"
                    # Label region with size threshold confidence
                    if pos2 - pos1 < int(REPORT_THRESHOLD):
                        confidence += "+low_confidence-size"
                    # Write to output bed file
                    outfile.write(f"{chrom}\t{pos1}\t{pos2}\t{confidence}\n")
                    # Add to combined DF
                    window_df_out = pd.concat([window_df_out, check_df.loc[y[0]:y[1]]])
                window_df_out.to_csv(output.windows_call, sep='\t', index=False)
