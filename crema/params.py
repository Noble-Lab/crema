"""
The :py:class:`Params` class is used to define the details and arguments necessary
for running crema from the command line.
"""

import argparse
import os


class Params:
    """
    All possible arguments and parameters for crema's CLI.
    Options can be specified as command-line arguments.
    """

    def __init__(self):
        """
        Initialize a Params object that holds an argparse parser.
        """
        self.parser = argparse.ArgumentParser(
            description="Calculate confidence estimates for PSMs"
        )
        self._configure_parser()

    def _configure_parser(self):
        """Configures all the arguments for the parser"""

        self.parser.add_argument(
            "input_files",
            type=str,
            nargs="+",
            help="One or more tab-delimited file(s) to read",
        )

        self.parser.add_argument(
            "--score",
            type=str,
            default="p-value",
            help="name of the column that defines the scores (p-values) of the psms. Defaults to 'p-value'",
        )

        self.parser.add_argument(
            "--spectrum",
            type=str,
            default="scan",
            help="name of the column that identifies the psm. Defaults to 'scan'",
        )

        self.parser.add_argument(
            "--target",
            type=str,
            default="target",
            help="name of the column that indicates if a psm is a target/decoy. Defaults to 'target'",
        )

        self.parser.add_argument(
            "--crux",
            action="store_true",
            help="Specifies that the input files are given in crux format",
        )

        self.parser.add_argument(
            "--file_root",
            type=str,
            help="This string will be added as a prefix to all output file names",
        )

        self.parser.add_argument(
            "--output_dir",
            type=str,
            default=os.getcwd(),
            help="The directory where output files will be created. Defaults to current working directory.",
        )

        self.parser.add_argument(
            "--logging",
            type=str,
            nargs="?",
            const=os.getcwd(),
            help="Specifies whether a logging file should be created. Defaults to current working directory if True.",
        )
