"""Writer to save confidence estimates in Parquet format."""

from pathlib import Path


def to_parquet(conf, output_dir=None, file_root=None, decoys=False):
    """Save confidence estimates to Parquet files.

    Writes one Parquet file per confidence level (PSMs, peptides, proteins,
    protein groups).  Requires the ``pyarrow`` package
    (``pip install crema[fast]``).

    Parameters
    ----------
    conf : Confidence object or tuple of Confidence objects
        One or more :py:class:`~crema.confidence.Confidence` objects.
    output_dir : str or None, optional
        Directory in which to save the files.  ``None`` uses the current
        working directory.
    file_root : str or None, optional
        Optional prefix for output file names.  Files are always named
        ``[file_root.]crema.{level}.parquet``.
    decoys : bool, optional
        Also write decoy confidence estimates.

    Returns
    -------
    list of str
        Paths to the saved files.
    """
    try:
        import pyarrow  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "The 'pyarrow' package is required to write Parquet files. "
            "Install it with: pip install crema[fast]"
        ) from exc

    from collections import defaultdict
    import pandas as pd
    from ..writers.txt import _get_level_data

    if isinstance(conf, str):
        raise ValueError("'conf' should be a Confidence object, not a string.")
    try:
        iter(conf)
    except TypeError:
        conf = [conf]

    file_base = "crema"
    if file_root is not None:
        file_base = file_root + "." + file_base
    if output_dir is not None:
        file_base = Path(output_dir, file_base)

    results = defaultdict(list)
    for res in conf:
        for level, qval_list in _get_level_data(res, decoys).items():
            results[level] += qval_list

    out_files = []
    for level, qval_list in results.items():
        out_file = str(file_base) + f".{level}.parquet"
        pd.concat(qval_list).to_parquet(out_file, index=False)
        out_files.append(out_file)

    return out_files
