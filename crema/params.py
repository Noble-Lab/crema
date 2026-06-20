"""The :py:class:`Params` class is used to define the details and arguments
necessary for running crema from the command line.
"""

import argparse
import textwrap

try:
    from . import __version__
except ImportError:
    __version__ = "unknown"


class CremaHelpFormatter(argparse.HelpFormatter):
    """Format help text to keep newlines and whitespace"""

    def _fill_text(self, text, width, indent):
        text_list = text.splitlines(keepends=True)
        return "\n".join(
            _process_line(line, width, indent) for line in text_list
        )


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

    subparsers = parser.add_subparsers(
        dest="command",
        metavar="COMMAND",
    )

    # ------------------------------------------------------------------
    # assign-confidence subcommand (original behaviour)
    # ------------------------------------------------------------------
    ac = subparsers.add_parser(
        "assign-confidence",
        help="Assign FDR confidence estimates to PSMs.",
        formatter_class=CremaHelpFormatter,
    )

    ac.add_argument(
        "psm_files",
        type=str,
        nargs="+",
        help=(
            "One or more collection of peptide-spectrum matches (PSMs) in the "
            "mzTab, Tide tab-delimited, MSGF+ tsv, MSAmanda csv, Morpheus txt, "
            "generic delimited text, or Parquet format."
        ),
    )

    ac.add_argument(
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

    ac.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=0.01,
        help=(
            "The FDR threshold for accepting discoveries. Default is 0.01. "
            "If 'q-value' is chosen, then 'accept' column is replaced "
            "with 'crema q-value'."
        ),
    )

    ac.add_argument(
        "-p",
        "--pep_fdr_type",
        type=str,
        default="psm-only",
        choices=["psm-only", "peptide-only", "psm-peptide"],
        help="The peptide-level FDR estimation method to use. "
        "Default is 'psm-only'",
    )

    ac.add_argument(
        "-r",
        "--prot_fdr_type",
        type=str,
        default="best",
        choices=["best", "combine"],
        help="The protein-level FDR estimation method to use. Default is 'best'",
    )

    ac.add_argument(
        "-f",
        "--file_root",
        type=str,
        help="This string will be added as a prefix to all output file names.",
    )

    ac.add_argument(
        "-o",
        "--output_dir",
        type=str,
        help=(
            "The directory where output files will be created. Defaults to "
            "the current working directory."
        ),
    )

    ac.add_argument(
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

    ac.add_argument(
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

    ac.add_argument(
        "-m",
        "--method",
        type=str,
        default="tdc",
        choices=["tdc"],
        help="The confidence estimation method to use.",
    )

    ac.add_argument(
        "--parquet",
        action="store_true",
        default=False,
        help="Write output as Parquet files instead of tab-delimited text.",
    )

    # ------------------------------------------------------------------
    # convert subcommand
    # ------------------------------------------------------------------
    cv = subparsers.add_parser(
        "convert",
        help="Convert PSM files to Parquet format for faster repeat analyses.",
        formatter_class=CremaHelpFormatter,
    )

    cv.add_argument(
        "psm_files",
        type=str,
        nargs="+",
        help=(
            "One or more PSM files to convert. Supported input formats: "
            "Tide tab-delimited, MSGF+ tsv, MSAmanda csv, Comet, MSFragger, "
            "generic delimited text, mzTab, pepXML."
        ),
    )

    cv.add_argument(
        "--target-column",
        type=str,
        required=True,
        help="Column indicating whether a PSM is a target or decoy.",
    )

    cv.add_argument(
        "--spectrum-columns",
        type=str,
        nargs="+",
        required=True,
        help="Column(s) that together uniquely identify a spectrum.",
    )

    cv.add_argument(
        "--score-columns",
        type=str,
        nargs="+",
        required=True,
        help="Column(s) containing PSM scores.",
    )

    cv.add_argument(
        "--peptide-column",
        type=str,
        required=True,
        help="Column containing peptide sequences.",
    )

    cv.add_argument(
        "--protein-column",
        type=str,
        required=True,
        help="Column containing protein identifiers.",
    )

    cv.add_argument(
        "--protein-delim",
        type=str,
        default=",",
        help="Delimiter separating multiple protein IDs. Default is ','.",
    )

    cv.add_argument(
        "-o",
        "--output",
        type=str,
        default=None,
        help=(
            "Output Parquet file path. Defaults to the first input file name "
            "with its extension replaced by '.parquet'."
        ),
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
