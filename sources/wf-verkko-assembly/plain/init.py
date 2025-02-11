#!/usr/bin/env python3

import argparse as argp
import logging
import logging.config as logconf
import os
import pathlib as pl
import subprocess as sp
import sys


# parent of the script's parent is one level above
# the repository location
DEFAULT_ROOT = pl.Path(__file__).resolve().parent.parent

DEFAULT_CONST_MODULE = pl.Path(__file__).parent
DEFAULT_CONST_MODULE = DEFAULT_CONST_MODULE.joinpath(
    "workflow", "rules", "commons", "10_constants.smk"
)
DEFAULT_CONST_MODULE.resolve(strict=True)


def create_execution_environment(repo_folder, project_folder, conda_env_name):
    """Create Conda environments if any Conda-like executable
    is found on $PATH

    Args:
        repo_folder (pathlib.Path): This repository checkout location
        project_folder (pathlib.Path): Project folder, assumed to be
        one above the repo folder
        conda_env_name: Name of Conda environment to create

    Raises:
        RuntimeError: In dev only mode, a Conda-like executable
        must be available

        subprocess.CalledProcessError: Propagated if Conda
        environment cannot be created

    Returns:
        None: placeholder
    """

    logger = logging.getLogger(__name__)
    check_executables = ["mamba", "conda"]
    use_executable = None
    for executable in check_executables:
        try:
            _ = sp.check_call(
                [executable, "--version"],
                shell=False,
                stdout=sp.DEVNULL,
                stderr=sp.DEVNULL,
            )
            use_executable = executable
            break
        except sp.CalledProcessError as spe:
            logger.warning(f"Executable {executable} not available: {spe}")
    if use_executable is None:
        logger.error(
            "No executable available to create execution (conda) environment."
            " If that is correct, please restart this script with the option"
            " '--no-env' to skip this step."
        )
        raise RuntimeError("No Conda/Mamba executable in $PATH")
    logger.debug(f"Found Conda executable {use_executable} - creating environment...")
    # select Conda env yaml file for specified environment
    std_path = repo_folder.joinpath("workflow", "envs", f"{conda_env_name}_env.yaml")
    logger.debug(f"Searching for Conda env file at location: {std_path}")
    yaml_file = std_path.resolve(strict=True)
    env_prefix = yaml_file.stem
    logger.debug(
        f"Creating environment with prefix '{env_prefix}/' underneath path: {project_folder}"
    )
    env_path = project_folder.joinpath(env_prefix)
    logger.debug("Setting up the Conda environment may take a while...")
    call_args = [
        use_executable,
        "env",
        "create",
        "--quiet",
        "--force",
        "-f",
        str(yaml_file),
        "-p",
        str(env_path),
    ]
    try:
        proc_out = sp.run(call_args, shell=False, capture_output=True, check=False)
        proc_out.check_returncode()  # check after to get stdout/stderr
    except sp.CalledProcessError as spe:
        logger.error(f"Could not create Snakemake execution environment: {spe}")
        logger.error(f"\n=== STDOUT ===\n{proc_out.stdout.decode('utf-8')}")
        logger.error(f"\n=== STDERR ===\n{proc_out.stderr.decode('utf-8')}")
        raise
    return None


def setup_logging(project_dir, debug_mode, dev_only):
    """Setup logging to stderr stream and file

    Args:
        project_dir (pathlib.Path): Project folder,
        assumed to be one above repo location
        debug_mode (boolean): log verbose / debug
        dev_only (boolean): create only Conda dev
        environment, do not create init.log file

    Returns:
        pathlib.Path or os.devnull: log file location
    """

    base_level = "DEBUG" if debug_mode else "WARNING"
    log_file_location = project_dir.joinpath("init.log")
    if dev_only:
        log_file_location = os.devnull

    log_config = {
        "version": 1,
        "root": {"handlers": ["stream", "file"], "level": "DEBUG"},
        "handlers": {
            "stream": {
                "formatter": "default",
                "class": "logging.StreamHandler",
                "level": base_level,
                "stream": sys.stderr,
            },
            "file": {
                "formatter": "default",
                "class": "logging.FileHandler",
                "level": "INFO",
                "filename": log_file_location,
            },
        },
        "formatters": {
            "default": {
                "format": (
                    "%(asctime)s : "
                    "%(levelname)s | "
                    "%(funcName)s | "
                    "ln:%(lineno)d >> "
                    "%(message)s"
                ),
                "datefmt": "%Y-%m-%d %H:%M:%S",
            },
        },
    }
    logconf.dictConfig(log_config)
    return log_file_location


def parse_command_line():
    """Create command line parser

    Returns:
        argparse.Namespace: command line options
    """
    parser = argp.ArgumentParser()
    parser.add_argument(
        "--debug",
        action="store_true",
        default=False,
        help="Print log messages to stderr.",
        dest="debug",
    )
    parser.add_argument(
        "--root-path",
        "--root",
        "-r",
        type=lambda x: pl.Path(x).resolve(strict=False),
        default=DEFAULT_ROOT,
        dest="root_path",
        help="Specify the top-level (root) path under which"
        " the Conda environment and the working directory"
        f" hierarchy shall be created. Default: {DEFAULT_ROOT}",
    )
    parser.add_argument(
        "--constants",
        "-c",
        type=lambda x: pl.Path(x).resolve(strict=True),
        default=DEFAULT_CONST_MODULE,
        dest="constants",
        help="Specify the path to the Snakemake 'constants' module"
        " that contains the information which paths to create"
        f" inside the working directory. Default: {DEFAULT_CONST_MODULE}",
    )
    parser.add_argument(
        "--dev-no-wd",
        "--no-wd",
        action="store_true",
        default=False,
        help="Do not create the working directory hierarchy.",
        dest="dev_no_wd",
    )
    parser.add_argument(
        "--dev-no-env",
        "--no-env",
        action="store_true",
        default=False,
        help="Do not create the Conda environment.",
        dest="dev_no_env",
    )
    parser.add_argument(
        "--conda-env",
        "--env",
        "-e",
        type=str,
        choices=["exec", "dev"],
        default="exec",
        dest="conda_env",
        help="Specify the Conda environment to create: 'exec' [default] or 'dev'",
    )

    args = parser.parse_args()
    return args


def create_wd_folders(project_dir, std_paths):
    """Create folder hierarchy starting
    at Snakemake's future working directory

    Args:
        project_dir (pathlib.Path): Project folder,
        assumed to be one above repo location
        std_paths (list of pathlib.Path): standard
        paths to be created in the working directory

    Returns:
        None: placeholder
    """

    logger = logging.getLogger(__name__)
    logger.info("Creating Snakemake working directory structure")
    wd_toplevel = project_dir.joinpath("wd")
    wd_toplevel.mkdir(exist_ok=True, parents=True)

    for sub_path in std_paths:
        full_path = wd_toplevel.joinpath(sub_path)
        logger.info(f"Creating path {full_path}")
        full_path.mkdir(exist_ok=True, parents=True)

    return None


def _extract_directory_paths(module_path):
    """
    Args:
        module_path (pathlib.Path): Path to Snakemake constants module
    """
    # import needed for eval() of Paths
    import pathlib

    logger = logging.getLogger(__name__)

    paths = dict()
    ignore = None
    extracting = False
    logger.debug("Evaluating content of constants module")
    with open(module_path, "r") as module:
        for line in module:
            if line.strip().startswith("#"):
                continue
            elif not line.strip() and extracting:
                # reached end of class definition
                break
            elif line.strip().startswith("class ConstDirectories:"):
                extracting = True
                continue
            elif line.strip().startswith("_no_init"):
                ignore = eval(line.strip().split("=")[-1])
            elif extracting:
                path_name = line.strip().split(":")[0].strip()
                path = eval(line.strip().split("=")[-1])
                paths[path_name] = path
            else:
                pass

    logger.debug(
        f"Extracted a total of {len(paths)} paths, {len(ignore)}"
        " of which to be ignored"
    )
    assert ignore is not None
    return paths, ignore


def load_constant_paths(module_path):
    """Load all default paths from the
    (default) constants module. The module
    is considered the single source of truth
    on standard paths required inside
    the working directory.

    Args:
        module_path (pathlib.Path): Path to Snakemake commons constants module
    """
    logger = logging.getLogger(__name__)
    logger.debug(f"Loading constants from module {module_path}")
    paths, ignore = _extract_directory_paths(module_path)
    paths_to_create = []
    for path_name, path in paths.items():
        if path_name not in ignore:
            paths_to_create.append(path)
    logger.debug(f"{len(paths_to_create)} remaining paths after filtering")
    return paths_to_create


def main():
    """Main function

    Returns:
        integer: explicit 0 on success
    """
    args = parse_command_line()
    dev_mode = args.conda_env == "dev"
    args.root_path.mkdir(parents=True, exist_ok=True)
    log_file_location = setup_logging(args.root_path, args.debug, dev_mode)
    logger = logging.getLogger(__name__)
    repo_location = pl.Path(__file__).parent
    logger.info(f"Repository location: {repo_location}")
    logger.info(f"Project directory: {args.root_path}")
    logger.info(f"Log file location: {log_file_location}")
    if not args.dev_no_env:
        create_execution_environment(repo_location, args.root_path, args.conda_env)
    else:
        logger.info("Skipping Conda environment creation.")
    if not args.dev_no_wd:
        paths = load_constant_paths(args.constants)
        create_wd_folders(args.root_path, paths)
    else:
        logger.info("Skipping creating working directory hierarchy.")

    return 0


if __name__ == "__main__":
    main()
