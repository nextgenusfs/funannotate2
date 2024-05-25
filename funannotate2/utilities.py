import os
import re
import textwrap
import shutil
from contextlib import contextmanager
import subprocess
import sys
import multiprocessing
import time
import sys
import queue
import signal
import uuid
from concurrent.futures import ThreadPoolExecutor
from threading import Lock
import requests
import errno
from urllib.request import urlopen
import socket
import random
import json
from .config import augustus_species, busco_taxonomy


def execute(cmd, cwd="."):
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
    with open(filename, "r") as f:
        data = json.load(f)
    return data


def download(url, name, wget=False):
    if wget:
        # download with wget
        cmd = ["wget", "-O", name, "--no-check-certificate", "-t", "2", "-c", url]
        subprocess.call(cmd)
    else:
        file_name = name
        try:
            u = urlopen(url)
            f = open(file_name, "wb")
            block_sz = 8192
            while True:
                buffer = u.read(block_sz)
                if not buffer:
                    break
                f.write(buffer)
            f.close()
        except socket.error as e:
            if e.errno != errno.ECONNRESET:
                raise
            pass


def lookup_taxonomy(name):
    """
    Fetch taxonomy dictionary given an organism species name

    Parameters
    ----------
    name : str
        name of organism to fetch, ie 'Aspergillus nidulans'

    Returns
    -------
    response : dict
        returns False if not a valid name else returns dictionary of taxonomic levels/names

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
                "superkingdom",
                "kingdom",
                "phylum",
                "class",
                "order",
                "family",
                "genus",
                "species",
            ]:
                if i in result:
                    data[i] = result[i]["name"]
            return data
    except requests.exceptions.RequestException as e:
        print("ERROR in taxonomy lookup: {}".format(e))
        return False


def pretty_taxonomy(
    obj, levels=["superkingdom", "kingdom", "phylum", "order", "class", "family"]
):
    # traverse taxonomy object and return list of taxonomy
    tax = []
    for t in levels:
        if t in obj:
            tax.append(obj[t])
    return tax


def choose_best_augustus_species(query_tax):
    return best_taxonomy(query_tax, augustus_species)


def choose_best_busco_species(query_tax):
    return best_taxonomy(query_tax, busco_taxonomy)


def best_taxonomy(query, ref):
    query_list = pretty_taxonomy(
        query,
        levels=[
            "superkingdom",
            "kingdom",
            "phylum",
            "order",
            "class",
            "family",
            "genus",
            "species",
        ],
    )
    ref_tax = {}
    for k, v in ref.items():
        simpletax = pretty_taxonomy(
            v,
            levels=[
                "superkingdom",
                "kingdom",
                "phylum",
                "order",
                "class",
                "family",
                "genus",
                "species",
            ],
        )
        res = [i for i, j in zip(query_list, simpletax) if i == j]
        if len(res) > 0:
            pos = simpletax.index(res[-1])
            ref_tax[k] = (len(res), len(simpletax[pos:]))
    # what is the max value
    iMax = max(ref_tax.items(), key=lambda x: x[1][0])[1][0]
    # now get a list of all of those keys
    best = {}
    for k, v in ref_tax.items():
        if v[0] == iMax:
            best[k] = v
    # finally we want the minimal value for ref tax
    iMin = min(best.items(), key=lambda x: x[1][1])[1][1]
    keepers = []
    for k, v in best.items():
        if v[1] == iMin:
            keepers.append(k)
    # here we can just take random
    if len(keepers) > 1:
        return random.choice(keepers)
    elif len(keepers) == 1:
        return keepers[0]
    else:
        return []


def which_path(file_name):
    for path in os.environ["PATH"].split(os.pathsep):
        full_path = os.path.join(path, file_name)
        if os.path.exists(full_path) and os.access(full_path, os.X_OK):
            return full_path
    return None


def human_readable_size(size, decimal_places=2):
    for unit in ["B", "KiB", "MiB", "GiB", "TiB", "PiB"]:
        if size < 1024.0 or unit == "PiB":
            break
        size /= 1024.0
    return f"{size:.{decimal_places}f} {unit}"


def create_tmpdir(outdir, base="predict"):
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


def checkfile(input):
    if input and os.path.isfile(input):
        filesize = os.stat(input).st_size
        if int(filesize) < 1:
            return False
        else:
            return True
    elif input and os.path.islink(input):
        return True
    else:
        return False


@contextmanager
def process_handle(handle, mode="w"):
    """
    A simple context manager for subprocess.run, yields and closes a file handle if a string is given, otherwise returns
    one of subprocess.PIPE, subprocess.DEVNULL, or NONE. If handle == "STDOUT", it will return subprocess.STDOUT for
    for redirecting stderr to stdout.

    :param handle: what type of return - True for PIPE, False for DEVNULL, None for None. A string will either return a
                   file handle to the path specified, or if handle=="STDOUT", subprocess.STDOUT (for redirecting stderr)
    :param mode: mode to open the file handle, one of "w", "r", "a"
    :return: PIPE, DEVNULL, STDOUT, NONE, or a file handle
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
    Runs a command using subprocess.run and directs output and stderr. If output or error are captured, they will be
    written to logfile.debug if the command succeeds, or logfile.error if the command fails.

    :param cmd: command to run
    :type cmd: list[str] or str
    :param cwd: string specifying the directory to run the command in
    :param logfile: logger for debug and error messages
    :param stdout: determines where to direct stdout. True - stdout will be captured; False - sent to /dev/null;
                   None - write to stdout; str - write to a file
    :type stdout: bool, str, or None
    :param stderr: determines where to direct stderr. True - stderr will be captured; False - sent to /dev/null;
                   None - write to stderr; str - write to a file, or if the str == "STDOUT", redirect to stdout
    :type stderr: bool, str, or None
    :param stdin: if reads input from file str
    :type stdin: str or None
    :param only_failed: only write stdout/stderr to the logfile if the command fails (returns a non-zero exit status).
                       This only effects stdout/stderr if it's being captured (output/error == True)
    :type only_failed: bool
    :param raise_not_exit: raise a subprocess.CalledProcessError if the subprocess fails, rather than sys.exit(1).
                           Default is False
    :type raise_not_exit: bool
    """

    def logoutput(logname, process_result):
        # stdout/stderr can be None or '', both evaluate to False, but '' != None
        if process_result.stdout:
            logname(process.stdout)
        if process_result.stderr:
            logname(process.stderr)

    logfile.debug(" ".join(cmd))
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
    # the command here is the function you want to run
    # the argList is a list of argument to pass [arg1, arg2, arg3]

    # simple progress indicator callback function
    def _progress_indicator(future):
        global lock, tasks_total, tasks_completed, tasks_failed
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
                f"  Progress: {tasks_completed}/{tasks_total} complete, {tasks_failed} failed, {tasks_total-tasks_completed} remaining        \r"
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
    global lock, tasks_total, tasks_completed, tasks_failed
    lock = Lock()
    tasks_total = len(argList)
    tasks_completed = 0
    tasks_failed = 0

    results = []
    executor = ThreadPoolExecutor(max_workers=cpus + 4)
    signal.signal(signal.SIGINT, lambda sig, frame: _exit_threads(executor))

    # submit jobs to executor
    for a in argList:
        results.append(executor.submit(func, *a))

    # add callback for progress
    if progress:
        for future in results:
            future.add_done_callback(_progress_indicator)

    return results


def runProcessJob(function, inputList, cpus=2):
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
    verbosity=1,
):
    """
    Args:
        table: a list of lists of strings (one row is one list, all rows should be the same length)
        alignments: a string of L and R, indicating the alignment for each row
        max_col_width: values longer than this will be wrapped
        col_separation: the number of spaces between columns
        indent: the number of spaces between the table and the left side of the terminal
        row_colour: a dictionary of row indices and their colour names
        sub_colour: a dictionary of values to colour names for which the text colour will be set
        row_extra_text: a dictionary of row indices and extra text to display after the row
        leading_newline: if True, the function will print a blank line above the table
        subsequent_indent: this string will be added to the start of wrapped text lines
        return_str: if True, this function will return a string of the table instead of printing it
        header_format: the formatting (colour, underline, etc) of the header line
        hide_header: if True, the header is not printed
        fixed_col_widths: a list to specify exact column widths (automatic if not used)
        left_align_header: if False, the header will follow the column alignments
        bottom_align_header: if False, the header will align to the top, like other rows
        verbosity: the table will only be logged if the logger verbosity is >= this value
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
