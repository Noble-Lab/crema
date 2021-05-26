"""
This is the command line interface for crema
"""
import os
import sys
import time
import logging

from .parsers.crux import read_crux
from .parsers.mztab import read_mztab
from .params import Params


def main():
    """The CLI entry point"""
    start_time = time.time()

    # Creates the parser for parse args and reads in command line arguments
    args = Params()

    # Set up logging files
    log_file = "crema.log.txt"
    if args.file_root is not None:
        log_file = args.file_root + "." + log_file
    if args.output_dir is None:
        args.output_dir = os.getcwd()

    # Configure logging
    logging.basicConfig(
        filename=os.path.join(args.output_dir, log_file),
        filemode="w+",
        level=logging.INFO,
        format="[%(levelname)s] %(message)s",
    )

    # Write logs to stderr as well
    logging.getLogger().addHandler(logging.StreamHandler())

    logging.info("crema")
    logging.info("Written by Donavan See and William E Fondrie in the ")
    logging.info(
        "Department of Genome Sciences at the University of Washington."
    )
    logging.info("Command issued:")
    logging.info("%s", " ".join(sys.argv))
    logging.info("")
    logging.info("Starting Analysis")
    logging.info("=================")

    # Create dataset object
    try:
        psms = read_crux(args.psm_files)
    except:  # TODO: Add specific exceptions
        psms = read_mztab(args.psm_files)

    conf = psms.assign_confidence(
        score_column=args.score,
        eval_fdr=args.eval_fdr,
        method=args.method,
    )

    # Write result to file
    logging.info("Writing results...")
    conf.to_txt(output_dir=args.output_dir, file_root=args.file_root)

    # Calculate how long the confidence estimation took
    end_time = time.time()
    total_time = end_time - start_time
    logging.info("==== DONE! =====")
    logging.info("Wall Time: %.2fs", total_time)


if __name__ == "__main__":
    main()
