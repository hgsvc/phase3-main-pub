
rule add_deepvariant_merge_label:
    input:
        vcf = WORKDIR_EVAL.joinpath(
            "results/regions", "{sample}",
            "{sample}.asmerr-hifi-onlyPRI.dv-wg.vcf.gz"
        ),
        sizes = DIR_PROC.joinpath(
            "asm_label_qc", "assembly_size",
            "{sample}.sizes.txt"
        )
    output:
        bed = DIR_RES.joinpath(
            "asm_label_qc", "merge_tables", "by-sample",
            "{sample}", "{sample}.deepvariant.mrg-labels.bed"
        )
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    run:
        import pandas as pd

        vcf = pd.read_csv(
            input.vcf, sep="\t", comment="#", header=None,
            names=["seq", "vcf_start", "name", "ref", "alt", "qual", "filter", "info", "format", "sample"],
            usecols=["seq", "vcf_start", "ref", "alt"]
        )

        # the DeepVariant error calls comprise the whole assembly, e.g., also the rDNA,
        # which is not the case for any of the other annotations. Hence, filter the VCF.
        known = set(pd.read_csv(input.sizes, sep="\t", header=None, names=["seq", "size"])["seq"].values)
        vcf = vcf.loc[vcf["seq"].isin(known), :].copy()
        assert vcf.shape[0] > 0

        def process_vcf_row(row):
            start = row.vcf_start - 1
            alt_lens = [len(alt) for alt in row.alt.split(",")]
            max_alt_len = max(alt_lens)
            ref_len = len(row.ref)
            if ref_len == max_alt_len:
                # SNV error
                end = start + ref_len
            elif ref_len < max_alt_len:
                # read alignments suggest missing bases
                end = start + ref_len + 1
            elif ref_len > max_alt_len:
                # read alignments suggest expansion
                end = start + ref_len
            else:
                raise ValueError(row)
            record_len = end - start
            assert record_len > 0
            record = row.seq, start, end, f"DEEPVR::{record_len}"
            return record

        norm_records = (vcf.apply(process_vcf_row, axis=1)).tolist()
        norm_records = pd.DataFrame.from_records(norm_records, columns=["seq", "start", "end", "label"])
        norm_records.sort_values(["seq", "start", "end"], inplace=True)
        norm_records.to_csv(output.bed, sep="\t", header=False, index=False)
    # END OF RUN BLOCK
