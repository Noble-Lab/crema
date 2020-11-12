"""
The :py:class:`Params` class is used to define the details and arguments necessary
for running crema from the command line.
"""

import argparse
import os
import crema


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


def _configure_parser():
    """Creates and configures all the arguments for the parser"""

    desc = (
        f"crema version {crema.__version__}\n\n"
        "Written by Donavan See (seed99@cs.washington.edu)\n\n"
        "Department of Genome Sciences at the University of Washington\n\n"
        "Official code website: <link tbd>\n\n"
        "More documentation and examples: <link tbd>"
    )

    parser = argparse.ArgumentParser(description=desc,)

    parser.add_argument(
        "input_files",
        type=str,
        nargs="+",
        help="One or more tab-delimited file(s) to read",
    )

    parser.add_argument(
        "--score",
        type=str,
        default="p-value",
        help="name of the column that defines the scores (p-values) of the psms. Defaults to 'p-value'",
    )

    parser.add_argument(
        "--spectrum",
        type=str,
        default="scan",
        help="name of the column that identifies the psm. Defaults to 'scan'",
    )

    parser.add_argument(
        "--target",
        type=str,
        default="target",
        help="name of the column that indicates if a psm is a target/decoy. Defaults to 'target'",
    )

    parser.add_argument(
        "--crux",
        action="store_true",
        help="Specifies that the input files are given in crux format",
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

    parser.add_argument(
        "--logging",
        type=str,
        nargs="?",
        const=os.getcwd(),
        help="Specifies whether a logging file should be created. Defaults to current working directory if True.",
    )

    return parser
