
rule ncbi_fcs_adaptor_screening:
    input:
        sif = NCBI_FCS_ADAPTOR_SIF,
        sh_script = NCBI_FCS_ADAPTOR_SCRIPT,
        fasta = rules.merge_and_tag_asm_units.output.mrg_fasta
    output:
        report = DIR_PROC.joinpath(
            "05-preprocess", "ncbi_fcs",
            "{sample}.adaptor.wd", "fcs_adaptor_report.txt"
        )
    benchmark:
        DIR_RSRC.joinpath(
            "05-preprocess", "ncbi_fcs",
            "{sample}.adaptor.rsrc")
    log:
        DIR_LOG.joinpath(
            "05-preprocess", "ncbi_fcs",
            "{sample}.adaptor.log")
    resources:
        mem_mb = lambda wildcards, attempt: 8192 * attempt,
        time_hrs = lambda wildcards, attempt: attempt
    params:
        out_dir = lambda wildcards, output: pathlib.Path(output.report).parent,
        taxonomy = NCBI_FCS_ADAPTOR_TAXONOMY
    shell:
        "mkdir -p {params.out_dir}"
            " && "
        "{input.sh_script} --fasta-input {input.fasta} --output-dir {params.out_dir} "
        "--{params.taxonomy} --container-engine singularity --image {input.sif} &> {log}"


rule ncbi_fcs_gx_contamination_screening:
    input:
        sif = NCBI_FCS_GX_SIF,
        db = NCBI_FCS_GX_DB_PATH,
        py_script = NCBI_FCS_GX_SCRIPT,
        fasta = rules.merge_and_tag_asm_units.output.mrg_fasta
    output:
        report = DIR_PROC.joinpath(
            "05-preprocess", "ncbi_fcs", "{sample}.gx-contam.wd",
            f"{{sample}}.asm-mrg-tag.{NCBI_FCS_GX_TAX_ID}.fcs_gx_report.txt"
        ),
    benchmark:
        DIR_RSRC.joinpath(
            "05-preprocess", "ncbi_fcs",
            f"{{sample}}.asm-mrg.{NCBI_FCS_GX_TAX_ID}.gx-contam.rsrc"
        )
    log:
        DIR_LOG.joinpath(
            "05-preprocess", "ncbi_fcs",
            f"{{sample}}.asm-mrg.{NCBI_FCS_GX_TAX_ID}.gx-contam.rsrc"
        )
    conda: DIR_ENVS.joinpath("biotools", "ncbi_fcs.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: int((384 + 192 * attempt) * 1024),
        time_hrs = lambda wildcards, attempt: 23 * attempt
    params:
        out_dir = lambda wildcards, output: pathlib.Path(output.report).parent,
        tax_id = NCBI_FCS_GX_TAX_ID,
        db_name = NCBI_FCS_GX_DB_NAME
    shell:
        "export FCS_DEFAULT_IMAGE={input.sif} ; "
        "python3 {input.py_script} screen genome --fasta {input.fasta} "
        "--out-dir {params.out_dir} --gx-db {input.db}/{params.db_name} "
        "--tax-id {params.tax_id} &> {log}"


rule normalize_merge_ncbi_fcs_adaptor_report:
    input:
        report = rules.ncbi_fcs_adaptor_screening.output.report,
        fasta = rules.merge_and_tag_asm_units.output.mrg_fasta
    output:
        report = DIR_RES.joinpath(
            "reports", "contamination", "{sample}.asm-mrg.fcs-report-adaptor.norm.tsv"
        ),
        stats = DIR_RES.joinpath(
            "reports", "contamination", "{sample}.asm-mrg.fcs-report-adaptor.stats.tsv"
        )
    conda: DIR_ENVS.joinpath("pyutils.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt,
        time_hrs = lambda wildcards, attempt: attempt
    params:
        script = find_script("normalize_merge_report")
    shell:
        "{params.script} --adaptor --report {input.report} "
            "--fasta {input.fasta} --table {output.report} "
            "--statistics {output.stats}"


rule normalize_merge_ncbi_fcs_contamination_report:
    input:
        report = rules.ncbi_fcs_gx_contamination_screening.output.report,
        fasta = rules.merge_and_tag_asm_units.output.mrg_fasta
    output:
        report = DIR_RES.joinpath(
            "reports", "contamination", "{sample}.asm-mrg.fcs-report-gxcontam.norm.tsv"
        ),
        stats = DIR_RES.joinpath(
            "reports", "contamination", "{sample}.asm-mrg.fcs-report-gxcontam.stats.tsv"
        )
    conda: DIR_ENVS.joinpath("pyutils.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt,
        time_hrs = lambda wildcards, attempt: attempt
    params:
        script = find_script("normalize_merge_report")
    shell:
        "{params.script} --contamination --report {input.report} "
            "--fasta {input.fasta} --table {output.report} "
            "--statistics {output.stats}"


rule run_all_ncbi_fcs_reports:
    input:
        rep_adapter = expand(
            rules.normalize_merge_ncbi_fcs_adaptor_report.output.report,
            sample=SAMPLES,
        ),
        rep_gxcontam = expand(
            rules.normalize_merge_ncbi_fcs_contamination_report.output.report,
            sample=SAMPLES
        )
