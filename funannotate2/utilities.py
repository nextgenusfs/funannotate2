import json
import multiprocessing
import os
import queue
import random
import re
import shutil
import signal
import socket
import subprocess
import sys
import textwrap
import time
import uuid
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from threading import Lock
from urllib.request import urlopen

import requests

from .config import augustus_species, busco_taxonomy

# disable insecure warning
requests.packages.urllib3.disable_warnings()


def merge_coordinates(intervals):
    """
    Merge overlapping intervals from a list of intervals.

    This function takes a list of intervals, each represented as a list or tuple of two integers [start, end],
    and merges all overlapping intervals. The resulting list contains non-overlapping intervals sorted
    by their starting points.

    Parameters:
    - intervals (List[Union[List[int], Tuple[int, int]]]): A list of intervals to be merged.

    Returns:
    - List[List[int]]: A list of merged intervals.

    Example:
    merge_coordinates([[1, 3], [2, 6], [8, 10], [15, 18]]) -> [[1, 6], [8, 10], [15, 18]]
    merge_coordinates([(1, 3), (2, 6), (8, 10), (15, 18)]) -> [[1, 6], [8, 10], [15, 18]]
    """
    if not intervals:
        return []
    # Sort the intervals based on the starting points
    intervals = sorted(intervals, key=lambda x: x[0])
    # Convert first interval to a list to make it mutable
    merged_intervals = [[intervals[0][0], intervals[0][1]]]
    for current in intervals[1:]:
        last_merged = merged_intervals[-1]
        # If the current interval overlaps with the last merged interval, merge them
        if current[0] <= last_merged[1]:
            last_merged[1] = max(last_merged[1], current[1])
        else:
            # Otherwise, add the current interval to the merged list
            merged_intervals.append([current[0], current[1]])
    return merged_intervals


def execute(cmd, cwd="."):
    """
    Execute a shell command and yield its output line by line.

    This function runs a specified command in the shell within a given working directory.
    It yields each line of the command's standard output as it becomes available. If the
    command exits with a non-zero status, a `subprocess.CalledProcessError` is raised.

    Parameters:
    - cmd (str): The command to be executed.
    - cwd (str, optional): The working directory where the command will be executed (default is the current directory).

    Yields:
    - str: Each line of output from the executed command.

    Raises:
    - subprocess.CalledProcessError: If the command exits with a non-zero status.
    """
    DEVNULL = open(os.devnull, "w")
    popen = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, universal_newlines=True, stderr=DEVNULL, cwd=cwd
    )
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)


def naming_slug(species, strain, lowercase=False):
    """
    Generate a slug for a species and strain combination.

    This function creates a slug by combining the species and strain names. The species name
    is formatted by replacing spaces with underscores. The slug can be converted to lowercase
    if specified. If a strain is provided, it is appended to the species with underscores.

    Parameters:
    - species (str): The name of the species.
    - strain (str): The name of the strain.
    - lowercase (bool, optional): Whether to convert the slug to lowercase (default is False).

    Returns:
    - str: The generated slug combining species and strain.
    """
    combo = species.replace(" ", "_")
    if lowercase:
        combo = combo.lower()
    else:
        combo = combo.capitalize()
    if strain:
        slug = f"{combo}_{strain.replace(' ', '')}"
    else:
        slug = combo
    return slug


def load_json(filename):
    """
    Load JSON data from a file.

    This function opens a specified JSON file, reads its contents, and returns the data
    as a dictionary.

    Parameters:
    - filename (str): The path to the JSON file.

    Returns:
    - dict: The JSON data loaded from the file.
    """
    with open(filename, "r") as f:
        data = json.load(f)
    return data


def get_odb_version(downloads_json_file):
    """
    Extract BUSCO ODB version from the downloads JSON.

    This function parses the provided downloads JSON to extract the BUSCO ODB version
    for each specified taxonomy. It returns a dictionary mapping taxonomies to their
    corresponding ODB versions.

    Parameters:
    - downloads_json_file (file): The JSON file data containing download information.

    Returns:
    - str: The latest ODB version.
    """
    odb_versions = set()
    du = load_json(downloads_json_file)
    for taxon, info in du["busco"].items():
        if "_" in info[1]:
            odb_versions.add(info[1].rsplit("_", 1)[1])
    return sorted(odb_versions, reverse=True)[0]


def download(url, name, wget=False, timeout=60, retries=3):
    """
    Download a file from a given URL with improved error handling and retries.

    This function tries HTTPS first, then falls back to FTP if needed.

    Parameters:
    - url (str): The URL of the file to download.
    - name (str): The name to save the downloaded file as.
    - wget (bool, optional): Flag to use `wget` for downloading (default is False).
    - timeout (int, optional): Timeout in seconds for the download (default is 60).
    - retries (int, optional): Number of retry attempts for failed downloads (default is 3).

    Returns:
    None
    """
    if wget:
        # download with wget
        cmd = [
            "wget",
            "-O",
            name,
            "--no-check-certificate",
            "-t",
            str(retries),
            "-c",
            url,
        ]
        subprocess.call(cmd)
        return

    file_name = name
    attempt = 0

    # Try HTTPS first
    while attempt < retries:
        try:
            # Use requests for HTTP/HTTPS
            with requests.get(url, stream=True, timeout=timeout, verify=False) as r:
                r.raise_for_status()  # Raise an exception for HTTP errors
                with open(file_name, "wb") as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        if chunk:  # filter out keep-alive chunks
                            f.write(chunk)
            # If we get here, download was successful
            return
        except requests.exceptions.RequestException as e:
            attempt += 1
            if attempt < retries:
                print(
                    f"HTTPS download attempt {attempt} failed: {str(e)}. Retrying in 5 seconds..."
                )
                time.sleep(5)  # Wait before retrying
            else:
                print(f"HTTPS download failed after {retries} attempts: {str(e)}")

                # If HTTPS fails and the URL is using HTTPS to access an FTP server, try direct FTP
                if "ftp." in url and url.startswith("https://"):
                    print(f"Trying FTP fallback for {url}")
                    ftp_url = url.replace("https://", "ftp://")
                    ftp_success = _download_ftp(ftp_url, file_name, timeout)
                    if ftp_success:
                        return

                # If we get here, both HTTPS and FTP failed
                raise Exception(f"Download failed: {str(e)}")


def _download_ftp(url, file_name, timeout=60):
    """
    Helper function to download a file using FTP protocol.

    Parameters:
    - url (str): The FTP URL to download from.
    - file_name (str): The name to save the downloaded file as.
    - timeout (int, optional): Timeout in seconds for the download (default is 60).

    Returns:
    - bool: True if download was successful, False otherwise.
    """
    try:
        # Set a socket timeout for FTP connections
        socket.setdefaulttimeout(timeout)

        # Create a temporary file first to avoid partial downloads
        temp_file = file_name + ".tmp"

        u = None
        f = None
        try:
            u = urlopen(url)
            f = open(temp_file, "wb")
            block_sz = 8192
            while True:
                buffer = u.read(block_sz)
                if not buffer:
                    break
                f.write(buffer)

            # Close files before renaming
            if f:
                f.close()
                f = None
            if u:
                u.close()
                u = None

            # Rename temp file to final file
            if os.path.exists(file_name):
                os.remove(file_name)
            os.rename(temp_file, file_name)

            print(f"FTP download successful: {url}")
            return True

        except Exception as e:
            print(f"FTP download failed: {str(e)}")
            # Clean up temp file if it exists
            if f:
                f.close()
            if u:
                try:
                    u.close()
                except Exception:
                    pass
            if os.path.exists(temp_file):
                os.remove(temp_file)
            return False

    except Exception as e:
        print(f"FTP connection failed: {str(e)}")
        return False


def validate_augustus_species(species_name):
    """
    Validate that the provided Augustus species is available in the config.

    Parameters:
    - species_name (str): The Augustus species name to validate

    Returns:
    - bool: True if valid, False otherwise
    """
    return species_name in augustus_species


def validate_busco_lineage(lineage_name):
    """
    Validate that the provided BUSCO lineage is available in the config.

    Parameters:
    - lineage_name (str): The BUSCO lineage name to validate

    Returns:
    - bool: True if valid, False otherwise
    """
    return lineage_name in busco_taxonomy


def lookup_taxonomy(name):
    """
    Fetch taxonomy information for a given organism species name.

    This function retrieves a taxonomy dictionary for a specified organism by querying
    an external taxonomy service. If the species name is valid, it returns a dictionary
    containing taxonomic levels and names. If the name is invalid or an error occurs,
    it returns False. EDIT: JGI changed superkingdom to domain, so map it.

    Parameters:
    - name (str): The name of the organism to fetch, e.g., 'Aspergillus nidulans'.

    Returns:
    - dict: A dictionary of taxonomic levels and names, or False if the name is invalid.
    """
    if " " in name:
        name = name.replace(" ", "_")
    url = "https://taxonomy.jgi.doe.gov/stax/name/"
    try:
        resp = requests.get(url + name, verify=False)
        if 500 <= resp.status_code <= 599:
            return resp.status_code
        result = resp.json()[name]
        if "error" in result:
            return False
        else:
            # reformat result to simpler version
            data = {}
            for i in [
                "domain",
                "kingdom",
                "phylum",
                "class",
                "order",
                "family",
                "genus",
                "species",
            ]:
                if i in result:
                    if i == "domain":
                        data["superkingdom"] = result[i]["name"]
                    else:
                        data[i] = result[i]["name"]
            return data
    except requests.exceptions.RequestException as e:
        print("ERROR in taxonomy lookup: {}".format(e))
        return False


def pretty_taxonomy(obj, levels):
    """
    Extract specified taxonomy levels from a taxonomy object.

    This function traverses a taxonomy object and extracts the values for the specified
    taxonomy levels. It returns a list of these values in the order provided by the levels
    parameter.

    Parameters:
    - obj (dict): The taxonomy object to extract levels from.
    - levels (list, optional): A list of taxonomy levels to extract. Defaults to
      ["superkingdom", "kingdom", "phylum", "order", "class", "family"].

    Returns:
    - list: A list of taxonomy levels extracted from the object.
    """
    # traverse taxonomy object and return list of taxonomy
    return [obj.get(lvl) for lvl in levels]


def choose_best_augustus_species(query_tax):
    return best_taxonomy(query_tax, augustus_species)


def choose_best_busco_species(query_tax):
    return best_taxonomy(query_tax, busco_taxonomy, exact=True)


def best_taxonomy(query, reference, exact=False):
    """
    Find the best matching taxonomy in a reference dictionary based on a query taxonomy.

    This function compares the query taxonomy with reference taxonomies and returns the key
    of the taxonomy that best matches the query. It evaluates matches based on the number of
    matching taxonomy levels and selects the one with the minimal difference in remaining levels.

    Parameters:
    - query (dict): The query taxonomy object to compare.
    - ref (dict): The reference taxonomy dictionary to search for matches.
    - exact (boolean): If True, return the most specific match that fully defines a taxonomic level.

    Returns:
    - str or list: The key of the best matching taxonomy in the reference dictionary, or an empty list if no match is found.
    """
    levels = [
        "superkingdom",
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
    ]

    def normalize(value):
        return (
            value.lower() if isinstance(value, str) else value
        )  # Ensure strings are lowercase

    def similarity_score(query, ref):
        return sum(
            1
            for level in levels
            if level in query
            and level in ref
            and normalize(query[level]) == normalize(ref[level])
        )

    # Special case: If query fully defines a taxonomic level, return it directly
    if exact:
        for level in reversed(levels):  # Start from deepest and move upward
            if level in query and any(
                normalize(query[level]) == normalize(reference[name].get(level, ""))
                for name in reference
            ):
                return normalize(query[level])  # Return first valid match found

    best_matches = []
    highest_score = -1
    for name, attributes in reference.items():
        score = similarity_score(query, attributes)
        if score > highest_score:
            highest_score = score
            best_matches = [name]
        elif score == highest_score:
            best_matches.append(name)
    return random.choice(best_matches) if best_matches else None


def which_path(file_name):
    """
    Find the full path of a file by searching through the directories in the system's PATH environment variable.

    This function iterates over each directory listed in the system's PATH environment variable,
    checking if the specified file exists and is executable. If found, it returns the full path
    to the file; otherwise, it returns None.

    Parameters:
    - file_name (str): The name of the file to search for.

    Returns:
    - str or None: The full path of the file if found and executable, otherwise None.
    """
    for path in os.environ["PATH"].split(os.pathsep):
        full_path = os.path.join(path, file_name)
        if os.path.exists(full_path) and os.access(full_path, os.X_OK):
            return full_path
    return None


def human_readable_size(size, decimal_places=2):
    """
    Convert a size in bytes to a human-readable string representation.

    This function takes a size in bytes and converts it to a more human-readable format,
    using units such as B, KiB, MiB, GiB, TiB, and PiB. The size is formatted to the specified
    number of decimal places.

    Parameters:
    - size (float): The size in bytes to convert.
    - decimal_places (int, optional): The number of decimal places to include in the output (default is 2).

    Returns:
    - str: A string representing the size in a human-readable format with the appropriate unit.
    """
    for unit in ["B", "KiB", "MiB", "GiB", "TiB", "PiB"]:
        if size < 1024.0 or unit == "PiB":
            break
        size /= 1024.0
    return f"{size:.{decimal_places}f} {unit}"


def create_tmpdir(outdir, base="predict"):
    """
    Create a temporary directory for storing files.

    This function generates a temporary directory with a unique name based on the provided
    base name and a UUID. If an output directory is specified, the temporary directory is
    created within it. If the output directory is "/tmp", a subdirectory is created; otherwise,
    the specified directory is used directly. If no output directory is provided, the temporary
    directory is created in the current working directory.

    Parameters:
    - outdir (str): The output directory path.
    - base (str, optional): The base name for the temporary directory (default is "predict").

    Returns:
    - str: The absolute path of the created temporary directory.
    """
    # create a tmpdir for some files
    tmpdirslug = "{}_{}".format(base, str(uuid.uuid4()))
    if outdir:
        if outdir == "/tmp":
            tmpdir = os.path.join(outdir, tmpdirslug)
        else:
            tmpdir = outdir
    else:
        tmpdir = tmpdirslug
    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)
    return os.path.abspath(tmpdir)


def create_directories(outdir, base="predict"):
    """
    Create necessary directories within the specified 'outdir' path.

    This function sets up a folder structure within the given output directory. If the 'outdir'
    does not exist, it creates the directory along with subdirectories for miscellaneous files,
    results, and log files. If 'outdir' already exists, it ensures the subdirectories are present
    and updates the 'base_results' directory if needed.

    Parameters:
    - outdir (str): The main directory path where the subdirectories will be created.
    - base (str, optional): The base name used for creating subdirectories (default is "predict").

    Returns:
    - tuple: A tuple containing the absolute paths to the created 'base_misc', 'base_results',
      and 'logfiles' directories.
    """
    # create folder structure
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
        os.makedirs(os.path.join(outdir, f"{base}_misc"))
        os.makedirs(os.path.join(outdir, f"{base}_results"))
        os.makedirs(os.path.join(outdir, "logfiles"))
    else:
        if os.path.isdir(os.path.join(outdir, f"{base}_results")):
            shutil.rmtree(os.path.join(outdir, f"{base}_results"))
            os.makedirs(os.path.join(outdir, f"{base}_results"))
        # make sure subdirectories exist
        dirs = [
            os.path.join(outdir, f"{base}_misc"),
            os.path.join(outdir, "logfiles"),
            os.path.join(outdir, f"{base}_results"),
        ]
        for d in dirs:
            if not os.path.isdir(d):
                os.makedirs(d)
    return (
        os.path.abspath(os.path.join(outdir, f"{base}_misc")),
        os.path.abspath(os.path.join(outdir, f"{base}_results")),
        os.path.abspath(os.path.join(outdir, "logfiles")),
    )


def find_files(directory, suffix):
    """
    Find files with a specific suffix in the given directory.

    This function searches through the specified directory and collects all files that
    have the given suffix. It returns a list of the full paths to these files.

    Parameters:
    - directory (str): The directory path to search for files.
    - suffix (str): The file suffix to filter files by.

    Returns:
    - list: A list of file paths that match the specified suffix in the directory.
    """
    hits = []
    for f in os.listdir(directory):
        if f.endswith(suffix):
            hits.append(os.path.join(directory, f))
    return hits


def checkfile(infile):
    """
    Check if the input file exists and is not empty.

    This function verifies whether the specified file exists and is not empty. It also checks
    if the input is a symbolic link. The function returns True if the file exists and is not
    empty or if it is a symbolic link; otherwise, it returns False.

    Parameters:
    - input (str): Path to the file to be checked.

    Returns:
    - bool: True if the file exists and is not empty or is a symbolic link, False otherwise.
    """
    if infile and os.path.isfile(infile):
        filesize = os.stat(infile).st_size
        if int(filesize) < 1:
            return False
        else:
            return True
    elif infile and os.path.islink(infile):
        return True
    else:
        return False


@contextmanager
def process_handle(handle, mode="w"):
    """
    Context manager for handling subprocess I/O streams.

    This function provides a context manager for managing file handles or standard I/O streams
    when using `subprocess.run`. It yields and closes a file handle if a string is provided,
    otherwise it returns one of `subprocess.PIPE`, `subprocess.DEVNULL`, or `None`. If the handle
    is "STDOUT", it returns `subprocess.STDOUT` for redirecting stderr to stdout.

    Parameters:
    - handle: Determines the type of return:
      - True for `subprocess.PIPE`
      - False for `subprocess.DEVNULL`
      - None for `None`
      - str for a file handle to the specified path, or `subprocess.STDOUT` if handle is "STDOUT".
    - mode (str): The mode to open the file handle, one of "w", "r", "a".

    Returns:
    - Yields `subprocess.PIPE`, `subprocess.DEVNULL`, `subprocess.STDOUT`, `None`, or a file handle.
    """
    close_on_exit = False  # Do we need to make sure the output is closed correctly?
    if handle is True:  # send to PIPE
        output = subprocess.PIPE
    elif handle is False:  # send to /dev/null
        output = subprocess.DEVNULL
    elif handle == "STDOUT":  # special case for piping stderr to stdout
        output = subprocess.STDOUT
    elif (
        handle is None
    ):  # use None (default for input/stdout/stderr) - generally, will write to stdout/stderr
        output = None
    else:
        mode_list = ["a", "r", "w"]
        if mode not in mode_list:
            raise ValueError(f"Mode must be one of {mode_list}, not {mode}")
        output = open(handle, mode)
        close_on_exit = True

    try:
        yield output
    finally:
        if close_on_exit:
            output.close()


def runSubprocess(
    cmd,
    logfile,
    cwd=".",
    stdout=True,
    stderr=True,
    stdin=None,
    only_failed=False,
    raise_not_exit=False,
    env=False,
):
    """
    Run a command using subprocess and direct output and stderr.

    This function executes a command using `subprocess.run` or `subprocess.Popen`, directing
    the standard output and error streams based on the provided parameters. If output or error
    streams are captured, they are logged to `logfile.debug` if the command succeeds, or
    `logfile.error` if it fails.

    Parameters:
    - cmd (list[str] or str): The command to run.
    - logfile: Logger for debug and error messages.
    - cwd (str, optional): Directory to run the command in (default is ".").
    - stdout (bool, str, or None, optional): Determines where to direct stdout:
      - True: Capture stdout.
      - False: Send to /dev/null.
      - None: Write to stdout.
      - str: Write to a file.
    - stderr (bool, str, or None, optional): Determines where to direct stderr:
      - True: Capture stderr.
      - False: Send to /dev/null.
      - None: Write to stderr.
      - str: Write to a file, or if "STDOUT", redirect to stdout.
    - stdin (str or None, optional): If provided, reads input from the specified file.
    - only_failed (bool, optional): Only write stdout/stderr to the logfile if the command fails (default is False).
    - raise_not_exit (bool, optional): Raise a `subprocess.CalledProcessError` if the subprocess fails, rather than exiting (default is False).
    - env (dict or None, optional): Environment variables to pass to the subprocess.

    Raises:
    - subprocess.CalledProcessError: If the command fails and `raise_not_exit` is True.
    - SystemExit: If the command fails and `raise_not_exit` is False.
    """

    def logoutput(logname, process_result):
        # stdout/stderr can be None or '', both evaluate to False, but '' != None
        if process_result.stdout:
            logname(process.stdout)
        if process_result.stderr:
            logname(process.stderr)

    try:
        logfile.debug(" ".join(cmd))
    except AttributeError:
        pass
    with process_handle(stdout) as p_out, process_handle(
        stderr
    ) as p_error, process_handle(stdin, mode="r") as p_in:
        if not env:
            process = subprocess.run(
                cmd,
                cwd=cwd,
                stdin=p_in,
                stdout=p_out,
                stderr=p_error,
                universal_newlines=True,
            )
        else:
            process = subprocess.Popen(
                cmd,
                cwd=cwd,
                stdin=p_in,
                stdout=p_out,
                stderr=p_error,
                universal_newlines=True,
                env=env,
            )
            process.communicate()
    try:
        process.check_returncode()  # Will raise a CalledProcessError if return code != 0
        if not only_failed:
            logoutput(logfile.debug, process)
    except subprocess.CalledProcessError as e:
        logfile.error(f"CMD ERROR: {' '.join(cmd)}")
        logoutput(logfile.error, process)
        if raise_not_exit:
            raise (e)
        else:
            raise SystemExit(1)


def runThreadJob(func, argList, cpus=2, progress=True):
    """
    Run multiple instances of a function concurrently using a thread pool.

    This function utilizes a thread pool to execute a given function concurrently with
    multiple sets of arguments. It provides an optional progress indicator to track the
    completion status of the tasks.

    Parameters:
    - func (callable): The function to run concurrently.
    - argList (list): A list of argument tuples to pass to the function.
    - cpus (int, optional): Number of CPUs to use (default is 2).
    - progress (bool, optional): Flag to show a progress indicator (default is True).

    Returns:
    - list: A list of futures representing the result of each function call.
    """
    # Global variables for thread job tracking
    lock = None
    tasks_total = 0

    # simple progress indicator callback function
    def _progress_indicator(future):
        global tasks_failed, tasks_completed
        # obtain the lock
        with lock:
            if future.cancelled():
                tasks_failed += 1
            elif future.exception():
                tasks_failed += 1
            else:
                # update the counter
                tasks_completed += 1
            # report progress
            sys.stdout.write(
                f"  Progress: {tasks_completed}/{tasks_total} complete, {tasks_failed} failed, {tasks_total - tasks_completed} remaining        \r"
            )
            if (
                tasks_total - tasks_completed == 0
                and tasks_completed + tasks_failed == tasks_total
            ):
                sys.stdout.write("\r")

    def _exit_threads(executor):
        print("\nWrapping up, please wait...")

        py_version = sys.version_info
        if (py_version.major == 3) and (py_version.minor < 9):
            # py versions less than 3.9
            # Executor#shutdown does not accept
            # cancel_futures keyword
            # manually shutdown
            # code taken from https://github.com/python/cpython/blob/3.9/Lib/concurrent/futures/thread.py#L210

            # prevent new tasks from being submitted
            executor.shutdown(wait=False)
            while True:
                # cancel all waiting tasks
                try:
                    work_item = executor._work_queue.get_nowait()

                except queue.Empty:
                    break

                if work_item is not None:
                    work_item.future.cancel()

        else:
            executor.shutdown(cancel_futures=True)

        sys.exit(0)

    # setup job here
    # create a lock for the counter
    lock = Lock()
    tasks_total = len(argList)
    tasks_completed = 0
    tasks_failed = 0

    results = []
    executor = ThreadPoolExecutor(max_workers=cpus + 4)
    signal.signal(signal.SIGINT, lambda _sig, _frame: _exit_threads(executor))

    # submit jobs to executor
    for a in argList:
        results.append(executor.submit(func, *a))

    # add callback for progress
    if progress:
        for future in results:
            future.add_done_callback(_progress_indicator)

    return results


def runProcessJob(function, inputList, cpus=2):
    """
    Run multiple instances of a function in parallel using multiprocessing.

    This function utilizes a multiprocessing pool to execute a given function in parallel
    across multiple CPUs. It processes a list of inputs, where each input is passed as
    arguments to the function. Results are collected and returned as a list.

    Parameters:
    - function (function): The function to be executed in parallel.
    - inputList (list): A list of lists, where each inner list represents an input to the function.
    - cpus (int, optional): The number of CPUs to be used for parallel processing (default is 2).

    Returns:
    - list: A list containing the results of each function execution.
    """

    # inputList here should be a list of lists with each input to function
    def update(res):
        results.append(res)

    def handle_error(error):
        print(f"Error: {error}", flush=True)

    # setup pool
    p = multiprocessing.Pool(cpus)
    # setup results and split over cpus
    results = []
    for i in inputList:
        if isinstance(i, list):
            p.apply_async(
                function, args=(i,), callback=update, error_callback=handle_error
            )
        else:
            p.apply_async(
                function, args=(i), callback=update, error_callback=handle_error
            )
    p.close()
    p.join()
    return results


def len_without_format(text):
    try:
        return len(remove_formatting(text))
    except TypeError:
        return len(str(text))


def remove_formatting(text):
    return re.sub("\033.*?m", "", text)


def colour(text, text_colour):
    bold_text = "bold" in text_colour
    text_colour = text_colour.replace("bold", "")
    underline_text = "underline" in text_colour
    text_colour = text_colour.replace("underline", "")
    text_colour = text_colour.replace("_", "")
    text_colour = text_colour.replace(" ", "")
    text_colour = text_colour.lower()
    if "red" in text_colour:
        coloured_text = RED
    elif "green" in text_colour:
        coloured_text = GREEN
    elif "yellow" in text_colour:
        coloured_text = YELLOW
    elif "dim" in text_colour:
        coloured_text = DIM
    else:
        coloured_text = ""
    if bold_text:
        coloured_text += BOLD
    if underline_text:
        coloured_text += UNDERLINE
    if not coloured_text:
        return text
    coloured_text += text + END_FORMATTING
    return coloured_text


END_FORMATTING = "\033[0m"
BOLD = "\033[1m"
UNDERLINE = "\033[4m"
RED = "\033[31m"
GREEN = "\033[32m"
MAGENTA = "\033[35m"
YELLOW = "\033[93m"
DIM = "\033[2m"


def green(text):
    return GREEN + text + END_FORMATTING


def bold_green(text):
    return GREEN + BOLD + text + END_FORMATTING


def red(text):
    return RED + text + END_FORMATTING


def magenta(text):
    return MAGENTA + text + END_FORMATTING


def bold_red(text):
    return RED + BOLD + text + END_FORMATTING


def bold(text):
    return BOLD + text + END_FORMATTING


def bold_underline(text):
    return BOLD + UNDERLINE + text + END_FORMATTING


def underline(text):
    return UNDERLINE + text + END_FORMATTING


def dim(text):
    return DIM + text + END_FORMATTING


def dim_underline(text):
    return DIM + UNDERLINE + text + END_FORMATTING


def bold_yellow(text):
    return YELLOW + BOLD + text + END_FORMATTING


def bold_yellow_underline(text):
    return YELLOW + BOLD + UNDERLINE + text + END_FORMATTING


def bold_red_underline(text):
    return RED + BOLD + UNDERLINE + text + END_FORMATTING


def print_table(
    table,
    alignments="",
    max_col_width=30,
    col_separation=3,
    indent=2,
    row_colour=None,
    sub_colour=None,
    row_extra_text=None,
    leading_newline=False,
    subsequent_indent="",
    return_str=False,
    header_format="underline",
    hide_header=False,
    fixed_col_widths=None,
    left_align_header=True,
    bottom_align_header=True,
):
    """
    Print a formatted table with customizable options.

    This function formats and prints a table based on the provided data and options. It supports
    alignment, column width customization, color formatting, and more. The table can be printed
    directly or returned as a string.

    Parameters:
    - table (list of list of str): A list of lists of strings, where each inner list represents a row.
    - alignments (str, optional): A string of 'L', 'R', or 'C' indicating the alignment for each column.
    - max_col_width (int, optional): Maximum width for each column; values longer than this will be wrapped (default is 30).
    - col_separation (int, optional): Number of spaces between columns (default is 3).
    - indent (int, optional): Number of spaces between the table and the left side of the terminal (default is 2).
    - row_colour (dict, optional): A dictionary of row indices and their color names.
    - sub_colour (dict, optional): A dictionary of values to color names for which the text color will be set.
    - row_extra_text (dict, optional): A dictionary of row indices and extra text to display after the row.
    - leading_newline (bool, optional): If True, the function will print a blank line above the table (default is False).
    - subsequent_indent (str, optional): String added to the start of wrapped text lines (default is "").
    - return_str (bool, optional): If True, the function will return a string of the table instead of printing it (default is False).
    - header_format (str, optional): The formatting (color, underline, etc.) of the header line (default is "underline").
    - hide_header (bool, optional): If True, the header is not printed (default is False).
    - fixed_col_widths (list, optional): A list to specify exact column widths (automatic if not used).
    - left_align_header (bool, optional): If False, the header will follow the column alignments (default is True).
    - bottom_align_header (bool, optional): If False, the header will align to the top, like other rows (default is True).
    - verbosity (int, optional): The table will only be logged if the logger verbosity is >= this value (default is 1).

    Returns:
    - str or None: Returns the formatted table as a string if `return_str` is True; otherwise, prints the table.
    """
    # this function is written by Ryan Wick in Unicycler code
    # modified to not support colors
    column_count = len(table[0])
    table = [x[:column_count] for x in table]
    table = [x + [""] * (column_count - len(x)) for x in table]
    if row_colour is None:
        row_colour = {}
    if sub_colour is None:
        sub_colour = {}
    if row_extra_text is None:
        row_extra_text = {}
    if leading_newline:
        print("")

    # Ensure the alignments string is the same length as the column count
    alignments += "L" * (column_count - len(alignments))
    alignments = alignments[:column_count]

    if fixed_col_widths is not None:
        col_widths = fixed_col_widths
    else:
        col_widths = [0] * column_count
        for row in table:
            col_widths = [
                min(max(col_widths[i], len_without_format(x)), max_col_width)
                for i, x in enumerate(row)
            ]
    separator = " " * col_separation
    indenter = " " * indent
    full_table_str = ""
    for i, row in enumerate(table):
        row = [str(x) for x in row]
        if hide_header and i == 0:
            continue

        if fixed_col_widths is not None:
            wrapped_row = []
            for col, fixed_width in zip(row, fixed_col_widths):
                wrapper = textwrap.TextWrapper(
                    subsequent_indent=subsequent_indent, width=fixed_width
                )
                wrapped_row.append(wrapper.wrap(col))
        else:
            wrapper = textwrap.TextWrapper(
                subsequent_indent=subsequent_indent, width=max_col_width
            )
            wrapped_row = [wrapper.wrap(x) for x in row]

        row_rows = max(len(x) for x in wrapped_row)
        if i == 0 and bottom_align_header:
            wrapped_row = [[""] * (row_rows - len(x)) + x for x in wrapped_row]

        for j in range(row_rows):
            row_line = [x[j] if j < len(x) else "" for x in wrapped_row]
            aligned_row = []
            for value, col_width, alignment in zip(row_line, col_widths, alignments):
                if alignment == "L" or (i == 0 and left_align_header):
                    aligned_row.append(value.ljust(col_width))
                elif alignment == "C":
                    aligned_row.append(value.center(col_width))
                else:
                    aligned_row.append(value.rjust(col_width))
            row_str = separator.join(aligned_row)
            if i in row_extra_text:
                row_str += row_extra_text[i]
            if i == 0 and header_format:
                row_str = colour(row_str, header_format)
            if i in row_colour:
                row_str = colour(row_str, row_colour[i])
            for text, colour_name in list(sub_colour.items()):
                row_str = row_str.replace(text, colour(text, colour_name))
            if j < row_rows - 1 and UNDERLINE in row_str:
                row_str = re.sub(r"\033\[4m", "", row_str)
            if return_str:
                full_table_str += indenter + row_str + "\n"
            else:
                print((indenter + row_str))
    if return_str:
        return full_table_str


def readBlocks(source, pattern):
    """
    Read lines from a source and yield blocks of lines that start with a specific pattern.

    This function processes a source line by line, grouping lines into blocks. A new block
    starts when a line begins with the specified `pattern`. Each block is yielded as a list
    of lines.

    Parameters:
    - source (iterable): The source to read lines from.
    - pattern (str): The pattern that the lines should start with to be considered as a new block.

    Yields:
    - list: A block of lines that start with the specified pattern.
    """
    buffer = []
    for line in source:
        try:
            line = line.decode("utf-8")
        except AttributeError:
            line = line
        if line.startswith(pattern):
            if buffer:
                yield buffer
            buffer = [line]
        else:
            buffer.append(line)
    yield buffer


def readBlocks2(source, startpattern, endpattern):
    """
    Read lines from a source and yield blocks of lines based on start and end patterns.

    This function processes a source line by line, grouping lines into blocks. A new block
    starts when a line begins with the specified `startpattern` or ends with the `endpattern`.
    Each block is yielded as a list of lines.

    Parameters:
    - source (iterable): The source to read line by line.
    - startpattern (str): The pattern that marks the start of a block.
    - endpattern (str): The pattern that marks the end of a block.

    Yields:
    - list: A block of lines that start with the `startpattern` or end with the `endpattern`.
    """
    buffer = []
    for line in source:
        try:
            line = line.decode("utf-8")
        except AttributeError:
            line = line
        if line.startswith(startpattern) or line.endswith(endpattern):
            if buffer:
                yield buffer
            buffer = [line]
        else:
            buffer.append(line)
    yield buffer


def rename_gff_contigs(gff, output, contigHeaderMap):
    """
    Rename contigs in a GFF3 file based on a provided mapping.

    This function reads a GFF3 file line by line and renames the contigs based on a provided
    mapping. The mapping is a dictionary where the keys are the original contig names and
    the values are the new names. The renamed GFF3 content is written to the specified output file.

    Parameters:
    - gff (str): The path to the input GFF3 file.
    - output (str): The path to the output GFF3 file.
    - contigHeaderMap (dict): A dictionary mapping original contig names to new names.

    Returns:
    - None
    """
    with open(output, "w") as outfile:
        with open(gff, "r") as infile:
            for line in infile:
                if line.startswith(">"):
                    continue
                if line.startswith("#"):
                    outfile.write(line)
                else:
                    cols = line.split("\t")
                    if cols[0] in contigHeaderMap:
                        cols[0] = contigHeaderMap[cols[0]]
                    outfile.write("\t".join(cols))
