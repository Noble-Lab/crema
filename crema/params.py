"""The :py:class:`Params` class is used to define the details and arguments
necessary for running crema from the command line.
"""
import argparse
from . import __version__


class CremaHelpFormatter(argparse.HelpFormatter):
    """Format help text to keep newlines and whitespace"""

    def _fill_text(self, text, width, indent):
        text_list = text.splitlines(keepends=True)
        return "\n".join(_process_line(l, width, indent) for l in text_list)


class Params:
    """
    All possible arguments and parameters for crema's CLI.
    Options can be specified as command-line arguments.
    """

    def __init__(self):
        """
        Initialize a Params object that holds an argparse parser.
        """
        self.parser = _configure_parser()
        self._namespace = vars(self.parser.parse_args())

    def __getattr__(self, option):
        return self._namespace[option]


def _configure_parser():
    """Creates and configures all the arguments for the parser"""

    desc = (
        f"crema version {__version__}\n\n"
        "Written by Donavan See (seed99@cs.washington.edu) and \n"
        "William E Fondrie (wfondrie@uw.edu) in the \n"
        "Department of Genome Sciences at the University of Washington\n\n"
        "Official code website: https://github.com/Noble-Lab/crema\n\n"
        "More documentation and examples: https://crema-ms.readthedocs.io/"
    )

    parser = argparse.ArgumentParser(
        description=desc, formatter_class=CremaHelpFormatter
    )

    parser.add_argument(
        "psm_files",
        type=str,
        nargs="+",
        help=(
            "One or more collection of peptide-spectrum matches (PSMs) in a "
            "tab-delimited format."
        ),
    )

    parser.add_argument(
        "--score",
        type=str,
        nargs="+",
        default="combined p-value",
        help="name of the column that defines the scores (p-values) of the psms."
        "\n Expects decimal or float column values",
    )

    parser.add_argument(
        "--spectrum",
        type=str,
        nargs="+",
        default="scan",
        help="one or more column names that identify the psms."
        "\n Expects numeric or string column values",
    )

    parser.add_argument(
        "--target",
        type=str,
        default="target/decoy",
        help="name of the column that indicates if a psm is a target/decoy."
        "\n Expects column values containing any of the following combinations: (True/False), (target/decoy),"
        "(t/d), (t/f), (1/0), (1/-1)",
    )

    parser.add_argument(
        "--score_choice",
        type=str,
        default="0",
        help="index of the score column to use for confidence estimation."
        "\n Expects column name as a string, or as an integer index from 0 to n:"
        "n being the total number of score columns.",
    )

    parser.add_argument(
        "--file_root",
        type=str,
        help="This string will be added as a prefix to all output file names",
    )

    parser.add_argument(
        "--output_dir",
        type=str,
        help="The directory where output files will be created. Defaults to current working directory.",
    )

    return parser


def _process_line(line, width, indent):
    """Process a line in the CLI help"""
    line = textwrap.fill(
        line,
        width,
        initial_indent=indent,
        subsequent_indent=indent,
        replace_whitespace=False,
    )
    return line.strip()
