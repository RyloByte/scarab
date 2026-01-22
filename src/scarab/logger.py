__author__ = 'Connor Morgan-Lang'

import logging
import os
import sys


class MyFormatter(logging.Formatter):
    """
    Controls the formatting of the log through the logging package
    """
    base_fmt = "%(asctime)s %(name)s %(levelname)s %(message)s"
    error_fmt = base_fmt
    warning_fmt = base_fmt
    debug_fmt = base_fmt
    info_fmt = base_fmt

    def __init__(self):
        super().__init__(fmt=MyFormatter.base_fmt,
                         datefmt="%d/%m %H:%M:%S")

    def format(self, record):

        # Save the original format configured by the user
        # when the logger formatter was instantiated
        format_orig = self._style._fmt

        # Replace the original format with one customized by logging level
        if record.levelno == logging.DEBUG:
            self._style._fmt = MyFormatter.debug_fmt

        elif record.levelno == logging.INFO:
            self._style._fmt = MyFormatter.info_fmt

        elif record.levelno == logging.ERROR:
            self._style._fmt = MyFormatter.error_fmt

        elif record.levelno == logging.WARNING:
            self._style._fmt = MyFormatter.warning_fmt

        # Call the original formatter class to do the grunt work
        result = logging.Formatter.format(self, record)

        # Restore the original format configured by the user
        self._style._fmt = format_orig
        if '\n' in result:
            message = record.getMessage()
            prefix = ''
            if message and message in result:
                prefix = result.split(message)[0]
            result = result.replace('\n', '\n' + prefix)

        return result


def prep_logging(log_file_name=None, verbosity=False):
    """
    Allows for multiple file handlers to be added to the root logger, but only a single stream handler.
    The new file handlers must be removed outside of this function explicitly
    :param log_file_name:
    :param verbosity:
    :return:
    """
    if verbosity:
        logging_level = logging.DEBUG
    else:
        logging_level = logging.INFO

    # Detect whether a handlers are already present and return if true
    logger = logging.getLogger()
    if len(logger.handlers):
        return

    formatter = MyFormatter()
    logger.setLevel(logging.DEBUG if log_file_name else logging_level)

    # Set the console handler normally writing to stdout/stderr
    ch = logging.StreamHandler()
    ch.setLevel(logging_level)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    if log_file_name:
        output_dir = os.path.dirname(log_file_name)
        try:
            if output_dir and not os.path.isdir(output_dir):
                os.makedirs(output_dir)
        except (IOError, OSError):
            sys.stderr.write("ERROR: Unable to make directory '" + output_dir + "'.\n")
            sys.exit(3)
        fh = logging.FileHandler(log_file_name, mode='w')
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)
        logger.addHandler(fh)
        logger.propagate = False
    return
