import sys
import os
import platform
import logging
import tracemalloc
from .__version__ import __version__
import numpy
import natsort
from .interlap import __version__ as interlap_version
from .utilities import human_readable_size
import gfftk
import buscolite


green = "\033[92m"
END = "\033[0m"


class CustomFormatter(logging.Formatter):
    green = "\033[92m"
    grey = "\x1b[38;20m"
    yellow = "\x1b[33;20m"
    red = "\x1b[31;20m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    end = "\033[0m"
    date_format = green + "%(asctime)s " + end
    msg_format = "%(message)s"
    file_format = "%(asctime)s %(message)s"

    FORMATS = {
        logging.DEBUG: date_format + grey + msg_format + reset,
        logging.INFO: date_format + grey + msg_format + reset,
        logging.WARNING: date_format + yellow + msg_format + reset,
        logging.ERROR: date_format + red + msg_format + reset,
        logging.CRITICAL: date_format + bold_red + msg_format + reset,
    }
    FILE_FORMATS = {
        logging.DEBUG: file_format,
        logging.INFO: file_format,
        logging.WARNING: file_format,
        logging.ERROR: file_format,
        logging.CRITICAL: file_format,
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt, datefmt="[%b %d %I:%M %p]")
        return formatter.format(record)


class FileFormatter(logging.Formatter):
    file_format = "%(asctime)s %(message)s"

    FILE_FORMATS = {
        logging.DEBUG: file_format,
        logging.INFO: file_format,
        logging.WARNING: file_format,
        logging.ERROR: file_format,
        logging.CRITICAL: file_format,
    }

    def format(self, record):
        log_fmt = self.FILE_FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt, datefmt="[%b %d %I:%M %p]")
        return formatter.format(record)


def startLogging(logfile=False):
    tracemalloc.start()
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    sth = logging.StreamHandler()
    sth.setLevel(logging.INFO)
    sth.setFormatter(CustomFormatter())
    logger.addHandler(sth)
    if logfile:
        # remove if exists else will just keep appending
        if os.path.isfile(logfile):
            os.remove(logfile)
        fhnd = logging.FileHandler(logfile)
        fhnd.setLevel(logging.DEBUG)
        fhnd.setFormatter(FileFormatter())
        logger.addHandler(fhnd)
    return logger


def system_info(log):
    log(
        "Python v{}; funannotate2 v{}; gfftk v{}; buscolite v{}".format(
            platform.python_version(),
            __version__,
            gfftk.__version__,
            buscolite.__version__,
        )
    )


def finishLogging(log, module):
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    log(
        "{} module finished: peak memory usage={}".format(
            module, human_readable_size(peak)
        )
    )


def log_dependencies(script=False):
    return 0
