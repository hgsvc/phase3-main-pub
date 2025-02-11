include: "rules/commons/00_commons.smk"
include: "rules/00_modules.smk"
include: "rules/99_aggregate.smk"

# To avoid adding a trivial testing
# module that extends the workflow
# output list, the WORKFLOW_OUTPUT
# for the snaketests is extended
# here. DO NOT DO this for a
# regular pipeline (use dedicated
# sub-modules).
WORKFLOW_OUTPUT.extend(
    [DIR_RES.joinpath("testing", "all-ok.txt")]
)

rule run_tests:
    input:
        RUN_CONFIG_RELPATH,
        COPY_SAMPLE_SHEET_RELPATH,
        MANIFEST_RELPATH,
        WORKFLOW_OUTPUT


rule run_tests_no_manifest:
    input:
        RUN_CONFIG_RELPATH,
        COPY_SAMPLE_SHEET_RELPATH,
        WORKFLOW_OUTPUT


rule create_test_file:
    """
    Implicitly tests pyutil functions
    02_pyutils.smk::get_username
    02_pyutils.smk::get_timestamp
    """
    output:
        DIR_PROC.joinpath("testing", "ts-ok_user-ok.txt"),
    params:
        acc_in=lambda wildcards, output: register_input(output, allow_non_existing=True),
    run:
        timestamp = get_timestamp()
        user_id = get_username()
        content = "Running Snakemake tests\n"
        content += f"get_username: {user_id}\n"
        content += f"get_timestamp: {timestamp}\n"
        with open(output[0], "w") as testfile:
            _ = testfile.write(content)
            _ = testfile.write("Creating test input file succeeded")
        # END OF RUN BLOCK



rule test_log_functions:
    """
    Test pyutil logging functions
    """
    output:
        DIR_PROC.joinpath("testing", "log-{logtype}-ok.txt"),
    params:
        acc_out=lambda wildcards, output: register_result(output),
    run:
        if wildcards.logtype == "out":
            logout("Test log message to STDOUT")
        elif wildcards.logtype == "err":
            logerr("Test log message to STDERR")
        else:
            raise ValueError(f"Unknown log type: {wildcards.logtype}")
        with open(output[0], "w") as testfile:
            testfile.write(f"Log test {wildcards.logtype} ok")
        # END OF RUN BLOCK



rule test_find_script_success:
    input:
        expand(rules.test_log_functions.output, logtype=["err", "out"]),
        rules.create_test_file.output,
    output:
        DIR_PROC.joinpath("testing", "success-find-script-ok.txt"),
    params:
        script=find_script("test"),
        acc_out=lambda wildcards, output: register_result(output),
    run:
        import pathlib  # workaround, see gh#20

        # the following should never raise,
        # i.e. script_find() would fail before
        _ = pathlib.Path(params.script).resolve(strict=True)
        with open(output[0], "w") as testfile:
            testfile.write("find_script success test ok")
        # END OF RUN BLOCK



rule test_find_script_fail:
    input:
        rules.test_find_script_success.output,
    output:
        DIR_PROC.joinpath("testing", "fail-find-script-ok.txt"),
    params:
        acc_out=lambda wildcards, output: register_result(output),
    run:
        try:
            script = find_script("non_existing")
            # the previous line should not succeed,
            # if we are here, we do not create the
            # output file of the test, and thus fail
        except ValueError:
            with open(output[0], "w") as testfile:
                testfile.write("find_script fail test ok")
        # END OF RUN BLOCK



rule test_rsync_f2d:
    input:
        rules.create_test_file.output,
        rules.test_find_script_fail.output,
    output:
        DIR_PROC.joinpath("testing", "subfolder", "ts-ok_user-ok.txt"),
    params:
        acc_out=lambda wildcards, output: register_result(output),
    run:
        import pathlib  # workaround, see gh#20

        # first check that nobody changed the filename
        input_name = pathlib.Path(input[0]).name
        output_name = pathlib.Path(output[0]).name
        assert input_name == output_name
        output_dir = pathlib.Path(output[0]).parent
        rsync_f2d(input[0], output_dir)
        # END OF RUN BLOCK



rule test_rsync_f2f:
    input:
        rules.create_test_file.output,
    output:
        DIR_PROC.joinpath("testing", "rsync-f2f-ok.txt"),
    params:
        acc_out=lambda wildcards, output: register_result(output),
    run:
        rsync_f2f(input[0], output[0])
        # END OF RUN BLOCK



rule test_rsync_fail:
    input:
        rules.create_test_file.output,
    output:
        DIR_PROC.joinpath("testing", "rsync-fail-ok.txt"),
    message:
        "EXPECTED FAILURE: ignore following rsync error message"
    params:
        acc_out=lambda wildcards, output: register_result(output),
    run:
        import subprocess

        try:
            rsync_f2d(input[0], "/")
        except subprocess.CalledProcessError:
            with open(output[0], "w") as testfile:
                testfile.write("rsync fail test ok")
        # END OF RUN BLOCK



rule test_git_labels:
    input:
        rules.test_rsync_f2d.output,
        rules.test_rsync_f2f.output,
        rules.test_rsync_fail.output,
    output:
        out=DIR_PROC.joinpath("testing", "git-labels-ok.txt"),
    params:
        acc_out=lambda wildcards, output: register_result(output.out),
    run:
        git_labels = collect_git_labels()
        with open(output[0], "w") as labels:
            for label, value in git_labels:
                _ = labels.write(f"{label}\t{value}\n")
        # END OF RUN BLOCK



if USE_REFERENCE_CONTAINER:
    CONTAINER_TEST_FILES = [
        DIR_GLOBAL_REF.joinpath("genome.fasta.fai"),
        DIR_GLOBAL_REF.joinpath("exclusions.bed"),
        DIR_GLOBAL_REF.joinpath("hg38_full.fasta.fai"),
        DIR_PROC.joinpath(".cache", "refcon", "refcon_manifests.cache"),
    ]
    REGISTER_REFERENCE_FILES = CONTAINER_TEST_FILES[:3]
else:
    CONTAINER_TEST_FILES = []
    REGISTER_REFERENCE_FILES = []


rule trigger_tests:
    input:
        rules.test_git_labels.output,
        CONTAINER_TEST_FILES,
    output:
        DIR_RES.joinpath("testing", "all-ok.txt"),
    params:
        acc_out=lambda wildcards, output: register_result(output),
        acc_ref=lambda wildcards, input: register_reference(REGISTER_REFERENCE_FILES),
    run:
        with open(output[0], "w") as testfile:
            testfile.write("ok")
        # END OF RUN BLOCK
