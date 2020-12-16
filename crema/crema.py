"""
This is the command line interface for crema
"""
import os
import sys
import time
import logging

from .parsers import *
from .params import Params
from .methods import calculate_tdc


def main():
    """The CLI entry point"""
    start_time = time.time()

    # Creates the parser for parse args and reads in command line arguments
    params = Params().parser
    args = params.parse_args()

    # Set up logging files
    log_file = "crema.logfile.log"
    if args.file_root is not None:
        log_file = args.file_root + log_file
    if args.output_dir is None:
        args.output_dir = os.getcwd()

    # Configure logging
    logging.basicConfig(
        filename=os.path.join(args.output_dir, log_file),
        level=logging.INFO,
        format="%(asctime)s %(message)s",
    )

    logging.info("crema")
    logging.info("Written by Donavan See (seed99@cs.washington.edu) in the")
    logging.info(
        "Department of Genome Sciences at the University of " "Washington."
    )
    logging.info("Command issued:")
    logging.info("%s", " ".join(sys.argv))
    logging.info("")
    logging.info("Starting Analysis")
    logging.info("=================")

    # Create dataset object
    logging.info("Creating dataset object...")
    psms = read_file(args.input_files, args.spectrum, args.score, args.target)

    # Run confidence estimate method
    logging.info("Calculating confidence estimate...")
    result = calculate_tdc(psms)

    # Write result to file
    logging.info("Writing to file...")
    result.write_file(output_dir=args.output_dir, file_root=args.file_root)

    # Calculate how long the confidence estimation took
    end_time = time.time()
    total_time = end_time - start_time

    logging.info("=== DONE! ===")
    logging.info("Time Taken:" + str(total_time))


if __name__ == "__main__":
    main()
