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
        "Official code website: https://github.com/Noble-Lab/crema\n\n"
        "More documentation and examples: https://crema-ms.readthedocs.io/"
    )

    parser = argparse.ArgumentParser(
        description=desc,
    )

    parser.add_argument(
        "input_files",
        type=str,
        nargs="+",
        help="One or more tab-delimited file(s) to read",
    )

    parser.add_argument(
        "--score",
        type=str,
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
