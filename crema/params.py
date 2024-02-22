"""The :py:class:`Params` class is used to define the details and arguments
necessary for running crema from the command line.
"""

import argparse
import textwrap
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
        "Written by Andy Lin (linandy@uw.edu), Donavan See (seed99@cs.washington.edu), "
        "and William E Fondrie (wfondrie@uw.edu) in the "
        "Department of Genome Sciences at the University of Washington "
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
            "One or more collection of peptide-spectrum matches (PSMs) in the "
            "mzTab, Tide tab-delimited formats, MSGF+ tsv file, MSAmanda csv "
            "output, Morpheus txt output, or generic delimited text format."
        ),
    )

    parser.add_argument(
        "-s",
        "--score",
        type=str,
        nargs="+",
        default=None,
        help=(
            "One or more columns that indicate possible scores by which to "
            "rank the PSMs. If more than one is provided, the best will be "
            "selected automatically. If none are provided, crema will try all "
            "available scores."
        ),
    )

    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=0.01,
        help=(
            "The FDR threshold for accepting discoveries. Default is 0.01. "
            "If 'q-value' is chosen, then “ accept”  column is replaced "
            "with 'crema q-value'."
        ),
    )

    parser.add_argument(
        "-p",
        "--pep_fdr_type",
        type=str,
        default="psm-only",
        choices=["psm-only", "peptide-only", "psm-peptide"],
        help="The peptide-level FDR estimation method to use."
        "Default is 'psm-peptide'",
    )

    parser.add_argument(
        "-r",
        "--prot_fdr_type",
        type=str,
        default="best",
        choices=["best", "combine"],
        help="The protein-level FDR estimation method to use. "
        "Default is 'best'",
    )

    parser.add_argument(
        "-f",
        "--file_root",
        type=str,
        help="This string will be added as a prefix to all output file names.",
    )

    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        help=(
            "The directory where output files will be created. Defaults to "
            "the current working directory."
        ),
    )

    parser.add_argument(
        "-d",
        "--desc",
        type=str,
        default="None",
        choices=["True", "False", "None"],
        help=(
            "True if higher scores better, False if lower scores are better. "
            "If None, crema will try both and use the choice that yields the "
            "most PSMs at the specified false discovery rate threshold "
            "(eval_fdr). If score_column is None, this parameter is ignored. "
            "Default is 'None'."
        ),
    )

    parser.add_argument(
        "-e",
        "--eval_fdr",
        type=float,
        default=0.01,
        help=(
            "The FDR threshold by which to choose the best score column and "
            "to report in logging messages. This should range from 0 to 1. "
            "Default value is 0.01."
        ),
    )

    parser.add_argument(
        "-m",
        "--method",
        type=str,
        default="tdc",
        choices=["tdc"],
        help="The confidence estimation method to use.",
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
