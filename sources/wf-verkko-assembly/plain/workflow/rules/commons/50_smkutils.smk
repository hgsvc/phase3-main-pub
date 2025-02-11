
localrules:
    dump_config,
    create_manifest,
    show_help,


rule show_help:
    """Rule to print a readable version
    of all supported config options.
    Dumps the same info into the file
    wd/help.generic-config.txt
    """
    output:
        "show-help",
    retries: 0
    message:
        "Dumping help info to file: help.generic-config.txt"
    run:
        import textwrap as twr
        import io
        import sys

        wrapper = twr.TextWrapper(
            width=80,
            initial_indent=4 * " ",
            subsequent_indent=6 * " ",
            fix_sentence_endings=True,
            break_on_hyphens=False,
        )
        # options can be omitted from help
        # to hide options irrelevant to
        # end users
        omit = None

        formatted_help = []
        for option, opt_spec in dataclasses.asdict(OPTIONS).items():
            if option == "omit":
                omit = set(opt_spec.default)
                continue
            wrapped_help = "\n".join(wrapper.wrap(opt_spec.help))
            opt_str = f" Option: {opt_spec.name}\n Help:\n{wrapped_help}\n\n"
            formatted_help.append(((option, opt_spec.name), opt_str))
        assert omit is not None

        buffer = io.StringIO()
        _ = buffer.write("\n====== WORKFLOW HELP =======")
        _ = buffer.write("\n== Generic config options ==\n\n")
        for (option, name), help_str in formatted_help:
            if option in omit or name in omit:
                continue
            _ = buffer.write(help_str)

        with open("help.generic-config.txt", "w") as dump:
            _ = dump.write(buffer.getvalue())

        sys.stdout.write(buffer.getvalue())
        # END OF RUN BLOCK



rule accounting_file_md5_size:
    """
    Compute MD5 checksum and file size
    for accounted files (inputs, references, results).
    For convenience, this rule also computes the file
    size to avoid this overhead at some other place
    of the pipeline.
    """
    input:
        source=load_file_by_path_id,
    output:
        md5=DIR_PROC.joinpath(
            ".accounting", "checksums", "{file_type}", "{file_name}.{path_id}.md5"
        ),
        file_size=DIR_PROC.joinpath(
            ".accounting", "file_sizes", "{file_type}", "{file_name}.{path_id}.bytes"
        ),
    benchmark:
        DIR_RSRC.joinpath(
            ".accounting", "checksums", "{file_type}", "{file_name}.{path_id}.md5.rsrc"
        )
    wildcard_constraints:
        file_type="(" + "|".join(sorted(ACCOUNTING_FILES.keys())) + ")",
    resources:
        time_hrs=lambda wildcards, attempt: 1 * attempt,
        mem_gb=lambda wildcards, attempt: 1 * attempt,
    shell:
        "md5sum {input.source} > {output.md5}"
        " && "
        "stat -c %s {input.source} > {output.file_size}"


rule accounting_file_sha256:
    """
    Compute SHA256 checksum (same as for MD5)
    """
    input:
        source=load_file_by_path_id,
    output:
        sha256=DIR_PROC.joinpath(
            ".accounting", "checksums", "{file_type}", "{file_name}.{path_id}.sha256"
        ),
    benchmark:
        DIR_RSRC.joinpath(
            ".accounting",
            "checksums",
            "{file_type}",
            "{file_name}.{path_id}.sha256.rsrc",
        )
    wildcard_constraints:
        file_type="(" + "|".join(sorted(ACCOUNTING_FILES.keys())) + ")",
    resources:
        time_hrs=lambda wildcards, attempt: 1 * attempt,
        mem_gb=lambda wildcards, attempt: 1 * attempt,
    shell:
        "sha256sum {input.source} > {output.sha256}"


rule dump_config:
    output:
        RUN_CONFIG_RELPATH,
    params:
        acc_in=lambda wildcards, output: register_input(output, allow_non_existing=True),
    run:
        import yaml

        runinfo = {"_timestamp": get_timestamp(), "_username": get_username()}
        git_labels = collect_git_labels()
        for label, value in git_labels:
            runinfo[f"_{label}"] = value
        # add complete Snakemake config
        runinfo.update(config)
        for special_key in ["devmode", "resetacc"]:
            try:
                del runinfo[special_key]
            except KeyError:
                pass

        with open(RUN_CONFIG_RELPATH, "w", encoding="ascii") as cfg_dump:
            yaml.dump(runinfo, cfg_dump, allow_unicode=False, encoding="ascii")
        # END OF RUN BLOCK



if SAMPLE_SHEET_NAME is not None:

    localrules:
        copy_sample_sheet,

    rule copy_sample_sheet:
        input:
            SAMPLE_SHEET_PATH,
        output:
            COPY_SAMPLE_SHEET_RELPATH,
        params:
            acc_in=lambda wildcards, output: register_input(
                output, allow_non_existing=True
            ),
        shell:
            "rsync {input} {output}"

else:

    localrules:
        no_sample_sheet,

    rule no_sample_sheet:
        """This is a mock-up rule
        needed because Snakemake cannot
        handle an empty input/output rule
        that would emerge above if no
        sample sheet is provided for the
        workflow run
        """
        output:
            COPY_SAMPLE_SHEET_RELPATH,
        shell:
            "touch {output}"


rule create_manifest:
    input:
        manifest_files=load_accounting_information,
    output:
        manifest=MANIFEST_RELPATH,
    run:
        import fileinput
        import collections
        import pandas

        # The following checks if accounting files are actually
        # in use - it's possible to write a workflow w/o using
        # reference files, and thus _not all_ of them have to
        # exist / be used. Part of fix for gh#15.
        process_accounting_files = {}
        for accounting_file, file_path in ACCOUNTING_FILES.items():
            if not file_path.is_file():
                if VERBOSE:
                    warn_msg = f"Warning: accounting file of type '{account_file}' not in use."
                    logerr(warn_msg)
                continue
            process_accounting_files[accounting_file] = file_path

        accounting_files_in_use = len(process_accounting_files)

        if accounting_files_in_use == 0:
            target_rule = "run_all_no_manifest"
            if NAME_SNAKEFILE == "snaketests":
                target_rule = "run_tests_no_manifest"

            err_msg = "No accounting files marked as in use.\n"
            err_msg += "This means one of three things:\n"
            err_msg += "0) You forgot to trigger the manifest creation\n"
            err_msg += "by running Snakemake in dry run mode twice\n"
            err_msg += "before the actual pipeline run.\n"
            err_msg += "1) Your workflow does not consume input, does not use\n"
            err_msg += "any reference file(s) and also does not produce output.\n"
            err_msg += "Really? Are you sure?\n"
            err_msg += "2) You did not annotate the workflow rules with:\n"
            err_msg += "commons/02_pyutils.smk::register_input()\n"
            err_msg += "commons/02_pyutils.smk::register_result()\n"
            err_msg += "commons/02_pyutils.smk::register_reference()\n"
            err_msg += "Please rerun the workflow twice in dry run mode...\n\n"
            err_msg += "snakemake --dry-run (or: -n) [...other options...]\n\n"
            err_msg += "...after fixing that.\n\n"
            err_msg += "However, if you are sure (!) that this is correct,\n"
            err_msg += f"please target the rule >>> {target_rule} <<<\n"
            err_msg += "to run the entire workflow w/o the manifest file.\n\n"
            logerr(err_msg)
            raise RuntimeError(
                "No accounts: cannot proceed with workflow execution w/o accouting files."
            )

        if len(input.manifest_files) == 0:
            assert accounting_files_in_use > 0
            # this combination of conditions can only
            # indicate an error
            err_msg = "No files collected to list in the MANIFEST, but\n"
            err_msg += f"{accounting_files_in_use} accounting files are\n"
            err_msg += "marked as in use.\n"
            err_msg += "Please check that you properly annotated rules\n"
            err_msg += "consuming input or reference files, and rules\n"
            err_msg += "producing output with the respective 'register_'\n"
            err_msg += "function from the commons/02_pyutils.smk module.\n\n"
            logerr(err_msg)
            raise RuntimeError(
                "No manifest files collected, but accounts are in use."
            )

        records = collections.defaultdict(dict)
        for line in fileinput.input(process_accounting_files.values(), mode="r"):
            path_id, path_record = process_accounting_record(line)
            records[path_id].update(path_record)

        df = pandas.DataFrame.from_records(list(records.values()))
        if df.empty:
            logerr("Manifest DataFrame is empty - aborting")
            raise RuntimeError("Manifest DataFrame is empty")

        df.sort_values(["file_category", "file_name"], ascending=True, inplace=True)
        reordered_columns = [
            "file_name",
            "file_category",
            "file_size_bytes",
            "file_checksum_md5",
            "file_checksum_sha256",
            "file_path",
            "path_id",
        ]
        assert all(c in df.columns for c in reordered_columns)
        df = df[reordered_columns]
        df.to_csv(output.manifest, header=True, index=False, sep="\t")
        # END OF RUN BLOCK
