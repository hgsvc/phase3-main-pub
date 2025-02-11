import datetime
import getpass
import pathlib
import pickle
import subprocess
import sys
import hashlib

# needed to handle reference container
# manifest cache files
import pandas

# locking capability potentially
# needed for file accounting;
# in the current implementation
# (only update during dry run)
# probably not needed
import portalocker


def logerr(msg):
    level = "ERROR"
    if VERBOSE:
        level = "VERBOSE (err/dbg)"
    write_log_message(sys.stderr, level, msg)
    return


def logout(msg):
    write_log_message(sys.stdout, "INFO", msg)
    return


def get_username():
    user = getpass.getuser()
    return user


def get_timestamp():
    # format: ISO 8601
    ts = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    return ts


def write_log_message(stream, level, message):
    # format: ISO 8601
    ts = get_timestamp()
    fmt_msg = f"{ts} - LOG {level}\n{message.strip()}\n"
    stream.write(fmt_msg)
    return


def find_script(script_name, extension="py"):

    predicate = lambda s: script_name == s.stem or script_name == s.name

    # DIR_SCRIPTS is set in common/constants
    all_scripts = DIR_SCRIPTS.glob(f"**/*.{extension.strip('.')}")
    retained_scripts = list(map(str, filter(predicate, all_scripts)))
    if len(retained_scripts) != 1:
        if len(retained_scripts) == 0:
            err_msg = (
                f"No scripts found or retained starting at '{DIR_SCRIPTS}' "
                f" and looking for '{script_name}' [+ .{extension}]"
            )
        else:
            ambig_scripts = "\n".join(retained_scripts)
            err_msg = f"Ambiguous script name '{script_name}':\n{ambig_scripts}\n"
        raise ValueError(err_msg)
    selected_script = retained_scripts[0]

    return selected_script


def flatten_nested_paths(struct):
    """
    Given an arbitrarily nested structure
    of items presumed to be Paths, return
    a flattened version containing only
    pathlib.Path objects.
    """

    flattened = set()
    for item in struct:
        if isinstance(item, str):
            flattened.add(pathlib.Path(item))
        elif isinstance(item, pathlib.Path):
            flattened.add(item)
        elif hasattr(item, "__iter__"):
            flattened = flattened.union(flatten_nested_paths(item))
        else:
            raise ValueError(f"Cannot handle item: {item}")
    return sorted(flattened)


# ==============================================
# Below: utility functions for file accounting
# ==============================================


def _add_file_to_account(accounting_file, fmt_records):
    """
    Add registered files to the respective accounting file;
    do not duplicate records by caching the set of known path
    ids in a separate (pickle dump) file.
    """
    lockfile = pathlib.Path(accounting_file).with_suffix(".lock")
    path_id_file = pathlib.Path(accounting_file).with_suffix(".paths.pck")
    with portalocker.Lock(lockfile, "a", timeout=WAIT_ACC_LOCK_SECS) as lock:
        try:
            with open(path_id_file, "rb") as path_id_dump:
                known_path_ids = pickle.load(path_id_dump)
        except (FileNotFoundError, EOFError):
            known_path_ids = set()

        with open(accounting_file, "a") as account:
            for path_id, path_records in fmt_records:
                if path_id in known_path_ids:
                    continue
                _ = account.write(path_records)
                known_path_ids.add(path_id)

        with open(path_id_file, "wb") as path_id_dump:
            _ = pickle.dump(known_path_ids, path_id_dump)
    return


def _add_checksum_files(add_files, subfolder):
    """
    Augment each file registered for manifest inclusion
    with three additional files recording md5 and sha256
    checksums and the file size in bytes.
    """
    formatted_records = []
    for add_file in add_files:
        file_name = add_file.name
        abspath = str(add_file.resolve(strict=False))
        path_id = hashlib.md5(abspath.encode("utf-8")).hexdigest()
        md5_checksum = DIR_PROC.joinpath(
            ".accounting", "checksums", f"{subfolder}", f"{file_name}.{path_id}.md5"
        )
        sha256_checksum = DIR_PROC.joinpath(
            ".accounting", "checksums", f"{subfolder}", f"{file_name}.{path_id}.sha256"
        )
        size_file = DIR_PROC.joinpath(
            ".accounting", "file_sizes", f"{subfolder}", f"{file_name}.{path_id}.bytes"
        )
        formatted_records.append(
            (
                path_id,
                f"{add_file}\t{path_id}\t{subfolder}\tdata\n"
                + f"{md5_checksum}\t{path_id}\t{subfolder}\tchecksum\n"
                + f"{sha256_checksum}\t{path_id}\t{subfolder}\tchecksum\n"
                + f"{size_file}\t{path_id}\t{subfolder}\tsize\n",
            )
        )

    return formatted_records


def register_input(*args, allow_non_existing=False):
    """
    This register function has slightly
    different semantics because input files
    should always exist when the workflow starts.
    For special cases, a keyword-only argument can
    be set to accept non-existing files when the
    pipeline run starts and yet the files should
    be counted as input files. One of those
    special cases is the config dump, which is
    counted as part of the input (cannot run
    the workflow w/o config), but the dump
    is only created after execution.

    Note: the return value must be constant to
    avoid that Snakemake triggers a rerun b/c
    of a changed rule parameter
    """
    if DRYRUN:
        accounting_file = ACCOUNTING_FILES["inputs"]
        input_files = flatten_nested_paths(args)
        if not allow_non_existing:
            err_msg = ""
            for file_path in input_files:
                if not file_path.exists():
                    err_msg += (
                        f"register_input() -> input file does not exist: {file_path}\n"
                    )
                elif file_path.is_dir():
                    err_msg += f"Specified input is not a regular file: {file_path}\n"
                else:
                    pass
            if err_msg:
                err_msg = f"\nINPUT ERROR:\n{err_msg}"
                logerr(err_msg)
                raise RuntimeError("Bad input files detected (see above)")
        fmt_records = _add_checksum_files(input_files, "inputs")
        _add_file_to_account(accounting_file, fmt_records)
    return None


def register_reference(*args):
    """
    See note in "register_input" regarding
    return value; same here.
    """
    if DRYRUN:
        accounting_file = ACCOUNTING_FILES["references"]
        reference_files = flatten_nested_paths(args)
        fmt_records = _add_checksum_files(reference_files, "references")
        _add_file_to_account(accounting_file, fmt_records)
    return None


def register_result(*args):
    """
    See note in "register_input" regarding
    return value; same here.
    """
    if DRYRUN:
        accounting_file = ACCOUNTING_FILES["results"]
        result_files = flatten_nested_paths(args)
        fmt_records = _add_checksum_files(result_files, "results")
        _add_file_to_account(accounting_file, fmt_records)
    return None


def load_accounting_information(wildcards):
    """
    This function loads the file paths from all
    three accounting files to force creation
    (relevant for checksum and size files).
    """
    created_files = []
    for account_name, account_file in ACCOUNTING_FILES.items():
        if VERBOSE:
            logerr(f"Loading file account {account_name}")
        try:
            with open(account_file, "r") as account:
                created_files.extend(
                    [l.split()[0] for l in account.readlines() if l.strip()]
                )
            if VERBOSE:
                logerr(f"Size of accounting file list: {len(created_files)}")
        except FileNotFoundError:
            if VERBOSE:
                warn_msg = f"Accounting file does not exist (yet): {account_file}\n"
                warn_msg += "Please RERUN the workflow in DRY RUN MODE to create the file accounts!"
                logerr(warn_msg)
    return sorted(created_files)


def load_file_by_path_id(wildcards):
    """ """
    account_type = wildcards.file_type
    assert account_type in ACCOUNTING_FILES
    account_file = ACCOUNTING_FILES[account_type]

    req_path_id = wildcards.path_id
    req_file = None
    with open(account_file, "r") as account:
        for line in account:
            # note here: the data file is always
            # first in the list of four entries
            file_path, file_path_id = line.split()[:2]
            if file_path_id != req_path_id:
                continue
            req_file = file_path
            break
    if req_file is None:
        logerr(
            f"Failed loading file with path ID {req_path_id} from accounting file {account_file}"
        )
        raise FileNotFoundError(
            f"Missing file: {wildcards.file_type} / {wildcards.file_name}"
        )

    return req_file


def _load_data_line(file_path):
    with open(file_path, "r") as dump:
        content = dump.readline().strip().split()[0]
    return content


def process_accounting_record(line):
    """ """
    path, path_id, file_category, record_type = line.strip().split()
    if record_type == "data":
        record = {
            "file_path": path,
            "file_name": pathlib.Path(path).name,
            "file_category": file_category,
            "path_id": path_id,
        }
    elif record_type == "checksum":
        checksum_type = path.rsplit(".", 1)[-1]
        record = {
            f"file_checksum_{checksum_type}": _load_data_line(path),
            "path_id": path_id,
        }
    elif record_type == "size":
        size_unit = path.rsplit(".", 1)[-1]
        record = {
            f"file_size_{size_unit}": int(_load_data_line(path)),
            "path_id": path_id,
        }
    return path_id, record


# ==============================================
# Below: utility functions for file copying
# ==============================================


def rsync_f2d(source_file, target_dir):
    """
    Convenience function to 'rsync' a source
    file into a target directory (file name
    not changed) in a 'run' block of a
    Snakemake rule.
    """
    abs_source = pathlib.Path(source_file).resolve(strict=True)
    abs_target = pathlib.Path(target_dir).resolve(strict=False)
    abs_target.mkdir(parents=True, exist_ok=True)
    _rsync(str(abs_source), str(abs_target))
    return


def rsync_f2f(source_file, target_file):
    """
    Convenience function to 'rsync' a source
    file to a target location (copy file
    and change name) in a 'run' block of
    a Snakemake rule.
    """
    abs_source = pathlib.Path(source_file).resolve(strict=True)
    abs_target = pathlib.Path(target_file).resolve(strict=False)
    abs_target.parent.mkdir(parents=True, exist_ok=True)
    _rsync(str(abs_source), str(abs_target))
    return


def _rsync(source, target):
    """
    Abstract function realizing 'rsync' calls;
    do not call this function, use 'rsync_f2f'
    or 'rsync_f2d'.
    """
    cmd = ["rsync", "--quiet", "--checksum", source, target]
    try:
        _ = subprocess.check_call(cmd, shell=False)
    except subprocess.CalledProcessError as spe:
        logerr(f"rsync from '{source}' to '{target}' failed")
        raise spe
    return


# ==============================================
# Below: utility functions for git interaction
# ==============================================


def _check_git_available():

    try:
        git_version = subprocess.check_output("git --version", shell=True)
        git_version = git_version.decode().strip()
        if VERBOSE:
            logerr(f"DEBUG: git version {git_version} detected in $PATH")
        git_in_path = True
    except subprocess.CalledProcessError:
        git_in_path = False
    return git_in_path


def _extract_git_remote():
    """
    Function to extract the (most likely) primary
    git remote name and URL (second part, split at ':').
    Applies following priority filter:

    1. github - CUBI convention
    2. githhu - CUBI convention
    3. origin - git default
    4. <other> (first in list, issues warning if VERBOSE is set)

    Returns:
        tuple of str: remote name, URL

    Raises:
        subprocess.CalledProcessError: if git executable
        is not available in PATH
        ValueError: if no git remotes are configured for the repo
    """

    remotes = []
    try:
        remotes = subprocess.check_output(
            "git remote -v", shell=True, cwd=DIR_REPOSITORY
        )
        remotes = remotes.decode().strip().split("\n")
        remotes = [tuple(r.split()) for r in remotes]
    except subprocess.CalledProcessError:
        error_msg = "ERROR:\n"
        error_msg = (
            f"Could not retrieve git remotes for repository: {DIR_REPOSITORY}\n"
            "Likely reason: this is not a git clone of the workflow repository.\n"
            "Consequence: the workflow repository information will not be recorded.\n"
            "This is strongly discouraged!"
        )
        # see gh#33 - no longer raises, although
        # this is really bad practice...
        logerr(error_msg)

    remote_name, remote_url = None, None

    if remotes:

        remote_priorities = {"github": 0, "githhu": 1, "origin": 2}

        # sort list of remotes by priority,
        # assign high rank / low priority to unexpected remotes
        remotes = sorted(
            [
                (remote_priorities.get(r[0], 10), r)
                for r in remotes
                if r[-1] == "(fetch)"
            ]
        )
        # drop priority value
        remote_name, remote_url, _ = remotes[0][1]
        remote_url = remote_url.split(":")[-1]

        if remote_name not in remote_priorities and VERBOSE:
            warning_msg = f"WARNING: unexpected git remote (name: {remote_name}) assumed to be primary."
            logerr(warning_msg)

    return remote_name, remote_url


def collect_git_labels():
    """
    Collect some basic information about the
    checked out git repository of the workflow
    """

    label_collection = {
        "git_remote": "unset-error",
        "git_url": "unset-error",
        "git_short": "unset-error",
        "git_long": "unset-error",
        "git_branch": "unset-error",
    }

    git_in_path = _check_git_available()

    if git_in_path:
        primary_remote, remote_url = _extract_git_remote()
    else:
        primary_remote, remote_url = None, None

    if primary_remote is not None:

        label_collection["git_remote"] = primary_remote
        label_collection["git_url"] = remote_url

        collect_options = [
            "rev-parse --short HEAD",
            "rev-parse HEAD",
            "rev-parse --abbrev-ref HEAD",
        ]
        info_labels = ["git_short", "git_long", "git_branch"]

        for option, label in zip(collect_options, info_labels):
            call = "git " + option
            try:
                # Important here to use DIR_SNAKEFILE (= the git repo location)
                # and not WORK_DIR, which would be the pipeline working directory.
                out = subprocess.check_output(call, shell=True, cwd=DIR_SNAKEFILE)
                out = out.decode().strip()
                assert label in label_collection
                label_collection[label] = out
            except subprocess.CalledProcessError as err:
                err_msg = f"\nERROR --- could not collect git info using call: {call}\n"
                err_msg += f"Error message: {str(err)}\n"
                err_msg += f"Call executed in path: {DIR_SNAKEFILE}\n"
                logerr(err_msg)

    git_labels = [(k, v) for k, v in label_collection.items()]

    return git_labels


# =======================================================
# Below: utility functions to handle reference container
# manifest caching
# =======================================================


def trigger_refcon_manifest_caching(wildcards):
    """
    This function merely triggers the checkpoint
    to merge all reference containers caches into
    one. This checkpoint is needed to get a
    start-to-end run, otherwise "refcon_find_container"
    would produce an error.
    """
    refcon_manifest_cache = str(
        checkpoints.refcon_cache_manifests.get(**wildcards).output.cache
    )
    expected_path = DIR_PROC.joinpath(".cache", "refcon", "refcon_manifests.cache")
    # following assert safeguard against future changes
    assert pathlib.Path(refcon_manifest_cache).resolve() == expected_path.resolve()
    return refcon_manifest_cache


def refcon_find_container(manifest_cache, ref_filename):

    if not pathlib.Path(manifest_cache).is_file():
        if DRYRUN:
            return ""
        else:
            if VERBOSE:
                warn_msg = "Warning: reference container manifest cache "
                warn_msg += "does not exist yet. Returning empty reference "
                warn_msg += "container path."
                logerr(warn_msg)
            return ""

    manifests = pandas.read_csv(manifest_cache, sep="\t", header=0)

    refcon_names = sorted(manifests["refcon_name"].unique())

    matched_names = set(manifests.loc[manifests["name"] == ref_filename, "refcon_name"])
    matched_alias1 = set(
        manifests.loc[manifests["alias1"] == ref_filename, "refcon_name"]
    )
    matched_alias2 = set(
        manifests.loc[manifests["alias2"] == ref_filename, "refcon_name"]
    )

    select_container = sorted(matched_names.union(matched_alias1, matched_alias2))
    if len(select_container) > 1:
        raise ValueError(
            f'The requested reference file name "{ref_filename}" exists in multiple containers: {select_container}'
        )
    elif len(select_container) == 0:
        raise ValueError(
            f'The requested reference file name "{ref_filename}" exists in none of these containers: {refcon_names}'
        )
    else:
        pass
    container_path = DIR_REFCON.joinpath(select_container[0] + ".sif")
    return container_path


def load_reference_container_names():

    existing_container = [sif_file.stem for sif_file in DIR_REFCON.glob("*.sif")]
    requested_container = config.get("reference_container_names", [])
    if not requested_container:
        raise ValueError(
            "The config option 'use_reference_container' is set to True. "
            "Consequently, you need to specify a list of container names "
            "in the config with the option 'reference_container_names'."
        )
    missing_container = ""
    for req_con in requested_container:
        if req_con not in existing_container:
            missing_container += f"\nMissing reference container: {req_con}"
            missing_container += f"\nExpected container image location: {DIR_REFCON.joinpath(req_con)}.sif\n"

    if missing_container:
        logerr(missing_container)
        raise ValueError(
            "At least one of the specified reference containers "
            "(option 'reference_container_names' in the config) "
            "does not exist in the reference container store."
        )
    return sorted(requested_container)


# ==============================================
# Below: utility functions for workflow staging
# ==============================================


def _extract_and_set_dryrun_constant():
    """
    The state of the '--dryrun' option is only available
    when parsing the original invocation command line.
    This function mimicks the respective code block
    in snakemake.__init__.py::main() (first few lines)
    """

    from snakemake import get_argument_parser as get_smk_cli_parser

    smk_cli_parser = get_smk_cli_parser()
    args, _ = smk_cli_parser.parse_known_args(sys.argv)

    cli_dryrun = args.dryrun
    assert isinstance(cli_dryrun, bool)

    # seems extremely unlikely to set the dryrun option as part
    # of an execution profile, so this just issues a warning if
    # VERBOSE is set
    if args.profile and VERBOSE:
        warning_msg = "\nWARNING (staging): the current value of option '--dryrun' "
        warning_msg += f"is set to {cli_dryrun}.\nThe '--profile' may override "
        warning_msg += "that value, but if so, this is ignored.\nWorkflow staging "
        warning_msg += f"proceeds with '--dry-run' set to {cli_dryrun}.\n"
        logerr(warning_msg)

    global DRYRUN
    DRYRUN = cli_dryrun

    return


def _reset_file_accounts():
    """
    In case erroneous entries in any of the
    file account prevent a proper execution
    of the pipeline, this function can be
    triggered by setting --config resetacc=True,
    and will delete the accounting files. Next,
    the workflow should be re-executed with
    `--dryrun` to recreate the accounting files.
    """
    if RESET_ACCOUNTING:
        for acc_name, acc_path in ACCOUNTING_FILES.items():
            if VERBOSE:
                logerr(f"Resetting accounting file {acc_name}")
            try:
                with open(acc_path, "w"):
                    pass
            except FileNotFoundError:
                pass
            # important: reset cached path IDs
            path_id_cache = pathlib.Path(acc_path).with_suffix(".paths.pck")
            try:
                with open(path_id_cache, "w"):
                    pass
            except FileNotFoundError:
                pass

    for acc_file_path in ACCOUNTING_FILES.values():
        acc_path = pathlib.Path(acc_file_path).parent
        acc_path.mkdir(exist_ok=True, parents=True)

    return
