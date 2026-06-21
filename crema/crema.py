"""
This is the command line interface for crema
"""

import os
import sys
import time
import logging
from pathlib import Path

from .parsers.tide import read_tide
from .parsers.msamanda import read_msamanda
from .parsers.msfragger import read_msfragger
from .parsers.msgf import read_msgf
from .parsers.comet import read_comet
from .parsers.mztab import read_mztab
from .parsers.pepxml import read_pepxml
from .parsers.txt import read_txt
from .params import Params


def main():
    """The CLI entry point"""
    start_time = time.time()

    args = Params()

    if args.command is None:
        args.parser.print_help()
        return

    if args.command == "assign-confidence":
        _run_assign_confidence(args, start_time)
    elif args.command == "convert":
        _run_convert(args, start_time)


def _setup_logging(output_dir, file_root):
    """Configure logging to file and stderr."""
    log_file = "crema.log.txt"
    if file_root is not None:
        log_file = file_root + "." + log_file
    if output_dir is None:
        output_dir = os.getcwd()

    logging.basicConfig(
        filename=os.path.join(output_dir, log_file),
        filemode="w+",
        level=logging.INFO,
        format="[%(levelname)s] %(message)s",
    )
    logging.getLogger().addHandler(logging.StreamHandler())

    logging.info("crema")
    logging.info(
        "Written by Andy Lin, Donavan See and William E Fondrie in the "
    )
    logging.info(
        "Department of Genome Sciences at the University of Washington."
    )
    logging.info("Command issued:")
    logging.info("%s", " ".join(sys.argv))
    logging.info("")
    logging.info("Starting Analysis")
    logging.info("=================")

    return output_dir


def _auto_read(psm_files):
    """Try each parser in turn; return the first that succeeds."""
    # Parquet is tried first via extension check to avoid slow fallback
    files = psm_files if isinstance(psm_files, list) else [psm_files]
    if all(str(f).endswith(".parquet") for f in files):
        raise ValueError(
            "Parquet input requires explicit column arguments. "
            "Use 'crema convert' output with 'crema assign-confidence'."
        )

    readers = [
        read_tide,
        read_msgf,
        read_msamanda,
        read_comet,
        read_msfragger,
        read_pepxml,
        read_mztab,
        read_txt,
    ]

    psms = None
    for read_fn in readers:
        try:
            psms = read_fn(psm_files)
            break
        except Exception as exc:
            logging.debug("%s failed: %s", read_fn.__name__, exc)
            continue

    if psms is None:
        raise ValueError("Unrecognized file type.")

    return psms


def _run_assign_confidence(args, start_time):
    """Run the assign-confidence subcommand."""
    output_dir = _setup_logging(
        getattr(args, "output_dir", None), getattr(args, "file_root", None)
    )

    # Check if inputs are Parquet
    files = args.psm_files
    if all(str(f).endswith(".parquet") for f in files):
        raise ValueError(
            "Parquet input to assign-confidence is not yet supported via the "
            "CLI. Load the Parquet file with crema.read_parquet() in Python."
        )

    psms = _auto_read(args.psm_files)

    desc = {"True": True, "False": False, "None": None}[args.desc]
    conf = psms.assign_confidence(
        score_column=args.score,
        threshold=args.threshold,
        pep_fdr_type=args.pep_fdr_type,
        prot_fdr_type=args.prot_fdr_type,
        desc=desc,
        eval_fdr=args.eval_fdr,
        method=args.method,
    )

    logging.info("Writing results...")
    if getattr(args, "parquet", False):
        conf.to_parquet(output_dir=output_dir, file_root=args.file_root)
    else:
        conf.to_txt(output_dir=output_dir, file_root=args.file_root)

    end_time = time.time()
    logging.info("==== DONE! =====")
    logging.info("Wall Time: %.2fs", end_time - start_time)


def _run_convert(args, start_time):
    """Run the convert subcommand."""
    logging.basicConfig(
        level=logging.INFO,
        format="[%(levelname)s] %(message)s",
    )

    logging.info("crema convert")
    logging.info("Command issued:")
    logging.info("%s", " ".join(sys.argv))
    logging.info("")

    logging.info("Reading PSMs...")
    psms = read_txt(
        args.psm_files,
        target_column=args.target_column,
        spectrum_columns=args.spectrum_columns,
        score_columns=args.score_columns,
        peptide_column=args.peptide_column,
        protein_column=args.protein_column,
        protein_delim=args.protein_delim,
    )

    if args.output is not None:
        out_path = args.output
    else:
        out_path = str(Path(args.psm_files[0]).with_suffix(".parquet"))

    logging.info("Writing Parquet to %s...", out_path)
    psms.data.to_parquet(out_path, index=False)

    end_time = time.time()
    logging.info("==== DONE! =====")
    logging.info("Wall Time: %.2fs", end_time - start_time)


if __name__ == "__main__":
    main()
