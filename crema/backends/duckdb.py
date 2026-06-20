"""DuckDB out-of-core backend for TDC confidence estimation.

This module implements the same target-decoy competition and q-value
calculation as :mod:`crema.qvalues`, but expressed as DuckDB SQL window
functions so that datasets larger than available RAM can be processed
without loading everything into memory at once.
"""

import logging

import pandas as pd

LOGGER = logging.getLogger(__name__)


def _require_duckdb():
    try:
        import duckdb

        return duckdb
    except ImportError as exc:
        raise ImportError(
            "The 'duckdb' package is required for backend='duckdb'. "
            "Install it with: pip install crema[large]"
        ) from exc


def _compete_into(con, src, dst, score_col, group_cols, desc):
    """Create table *dst* holding the best-scoring row per group from *src*.

    Ties are broken randomly using the ``_crema_rand`` column that must
    already exist in *src*.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
    src : str
        Name of an existing DuckDB table.
    dst : str
        Name of the table to create.
    score_col : str
        Column to rank by within each group.
    group_cols : list of str
        Columns that define a group (e.g. spectrum identifiers).
    desc : bool
        True if higher scores are better.
    """
    direction = "DESC" if desc else "ASC"
    partition = ", ".join(f'"{c}"' for c in group_cols)
    con.execute(
        f"""
        CREATE TABLE {dst} AS
        SELECT * EXCLUDE (_crema_rn)
        FROM (
            SELECT *,
                ROW_NUMBER() OVER (
                    PARTITION BY {partition}
                    ORDER BY "{score_col}" {direction}, _crema_rand
                ) AS _crema_rn
            FROM {src}
        )
        WHERE _crema_rn = 1
        """
    )


def _fetch_qvalues(con, table, score_col, target_col, desc):
    """Return a DataFrame from *table* with a ``crema q-value`` column added.

    Implements the same tied-score-aware, backward-monotone-minimum
    algorithm as :func:`crema.qvalues.tdc` using DuckDB window functions.

    The returned DataFrame is ordered worst-to-best (ascending score for
    ``desc=True``), matching the ordering that :class:`~crema.confidence.Confidence`
    expects before its ``_prettify_tables`` step reverses it.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
    table : str
        Name of the DuckDB table to read.
    score_col : str
    target_col : str
    desc : bool

    Returns
    -------
    pandas.DataFrame
    """
    # cum_dir: direction for cumulative sums (best → worst)
    # qval_dir: direction for running-min q-value pass (worst → best)
    cum_dir = "DESC" if desc else "ASC"
    qval_dir = "ASC" if desc else "DESC"

    return con.execute(
        f"""
        WITH cum_counts AS (
            SELECT *,
                SUM(CAST("{target_col}" AS INT)) OVER (
                    ORDER BY "{score_col}" {cum_dir}
                    ROWS UNBOUNDED PRECEDING
                ) AS _cum_t,
                SUM(CAST(NOT "{target_col}" AS INT)) OVER (
                    ORDER BY "{score_col}" {cum_dir}
                    ROWS UNBOUNDED PRECEDING
                ) AS _cum_d
            FROM {table}
        ),
        group_fdr AS (
            SELECT *,
                CASE
                    WHEN MAX(_cum_t) OVER (PARTITION BY "{score_col}") = 0
                        THEN 1.0
                    ELSE (MAX(_cum_d) OVER (PARTITION BY "{score_col}") + 1.0)
                         / MAX(_cum_t) OVER (PARTITION BY "{score_col}")
                END AS _group_fdr
            FROM cum_counts
        )
        SELECT * EXCLUDE (_cum_t, _cum_d, _group_fdr),
            MIN(_group_fdr) OVER (
                ORDER BY "{score_col}" {qval_dir}
                ROWS UNBOUNDED PRECEDING
            ) AS "crema q-value"
        FROM group_fdr
        ORDER BY "{score_col}" {qval_dir}
        """
    ).df()


def run_tdc(psms, score_column, desc, pep_fdr_type, prot_fdr_type, eval_fdr):
    """Run the full TDC pipeline using DuckDB as the compute engine.

    Replicates the logic of :meth:`crema.confidence.TdcConfidence._assign_confidence`
    using SQL window functions. Protein-group estimation is not implemented
    in this backend; use :class:`~crema.confidence.TdcConfidence` if protein
    groups are required.

    Parameters
    ----------
    psms : crema.dataset.PsmDataset
        The dataset to process.
    score_column : str
        Column to rank PSMs by.
    desc : bool
        True if higher scores are better.
    pep_fdr_type : {"psm-only", "peptide-only", "psm-peptide"}
        Method for peptide-level FDR.
    prot_fdr_type : {"best", "combine"}
        Method for protein score aggregation.
    eval_fdr : float
        FDR threshold used only for logging discovery counts.

    Returns
    -------
    tuple of (dict, dict)
        ``(confidence_estimates, decoy_confidence_estimates)`` — each maps
        level names to :class:`pandas.DataFrame` objects carrying a
        ``"crema q-value"`` column, split into targets and decoys.
    """
    duckdb = _require_duckdb()

    pairing = psms.peptide_pairing
    if pairing is None and pep_fdr_type != "psm-only":
        raise ValueError(
            "Must provide paired target decoy peptide information (see FAQ)."
        )

    target_col = psms._target_column
    spectrum_cols = psms._spectrum_columns
    pep_col = psms._peptide_column
    prot_col = psms._protein_column
    prot_delim = psms._protein_delim

    con = duckdb.connect()

    # Register the full dataset and materialise with a per-row random value
    # used for tie-breaking inside window ORDER BY clauses.
    data = psms.data.copy()
    con.register("_raw_data", data)
    con.execute(
        "CREATE TABLE raw_psms AS "
        "SELECT *, random() AS _crema_rand FROM _raw_data"
    )

    conf = {}
    decoy_conf = {}

    # ------------------------------------------------------------------
    # PSM level — best PSM per spectrum
    # ------------------------------------------------------------------
    LOGGER.info("DuckDB: computing PSM-level confidence...")
    LOGGER.warning(
        "PSM-level FDR estimates are not guaranteed to control the FDR. "
        "We suggest using peptide-level FDR estimates (see FAQ)."
    )
    _compete_into(
        con, "raw_psms", "competed_psms", score_column, spectrum_cols, desc
    )
    psm_df = _fetch_qvalues(con, "competed_psms", score_column, target_col, desc)
    psm_df = psm_df.drop(columns=["_crema_rand"])
    mask = psm_df[target_col].values.astype(bool)
    conf["psms"] = psm_df[mask].reset_index(drop=True)
    decoy_conf["psms"] = psm_df[~mask].reset_index(drop=True)
    LOGGER.info(
        "  - Found %i PSMs at q<=%g.",
        (conf["psms"]["crema q-value"] <= eval_fdr).sum(),
        eval_fdr,
    )

    # ------------------------------------------------------------------
    # Peptide level
    # ------------------------------------------------------------------
    LOGGER.info("DuckDB: computing peptide-level confidence...")
    if pep_fdr_type == "psm-only":
        # Compete over peptide column across all PSMs
        _compete_into(
            con, "raw_psms", "competed_pep", score_column, [pep_col], desc
        )
    else:
        # psm-peptide: start from spectrum-competed PSMs
        # peptide-only: start from all PSMs
        pep_src = "competed_psms" if pep_fdr_type == "psm-peptide" else "raw_psms"

        # Load the pairing map (target_seq → decoy_seq)
        pairing_df = pd.DataFrame(
            list(pairing.items()), columns=["_seq", "_pair"]
        )
        con.register("_pairing_data", pairing_df)
        con.execute(
            "CREATE TABLE _pairing AS SELECT * FROM _pairing_data"
        )

        # Replace the peptide column with its paired counterpart so that
        # target/decoy pairs compete together.  Decoys (not in the map) keep
        # their own sequence.
        con.execute(
            f"""
            CREATE TABLE with_pairing AS
            SELECT p.* EXCLUDE ("{pep_col}"),
                COALESCE(m._pair, p."{pep_col}") AS "{pep_col}"
            FROM {pep_src} p
            LEFT JOIN _pairing m ON p."{pep_col}" = m._seq
            """
        )
        _compete_into(
            con, "with_pairing", "competed_pep", score_column, [pep_col], desc
        )

    pep_df = _fetch_qvalues(con, "competed_pep", score_column, target_col, desc)
    pep_df = pep_df.drop(columns=["_crema_rand"])
    mask_pep = pep_df[target_col].values.astype(bool)
    conf["peptides"] = pep_df[mask_pep].reset_index(drop=True)
    decoy_conf["peptides"] = pep_df[~mask_pep].reset_index(drop=True)
    LOGGER.info(
        "  - Found %i Peptides at q<=%g.",
        (conf["peptides"]["crema q-value"] <= eval_fdr).sum(),
        eval_fdr,
    )

    # ------------------------------------------------------------------
    # Protein level — filter ambiguous proteins, aggregate, compete
    # ------------------------------------------------------------------
    LOGGER.info("DuckDB: computing protein-level confidence...")

    # Escape any single quotes in the delimiter for the SQL literal.
    safe_delim = prot_delim.replace("'", "''")
    con.execute(
        f"""
        CREATE TABLE unambig_psms AS
        SELECT * FROM competed_psms
        WHERE NOT CONTAINS("{prot_col}", '{safe_delim}')
        """
    )

    if prot_fdr_type == "best":
        agg_expr = f'{"MAX" if desc else "MIN"}("{score_column}")'
    else:
        # combine: sum of scores (desc) or product via log trick (asc)
        if desc:
            agg_expr = f'SUM("{score_column}")'
        else:
            agg_expr = f'EXP(SUM(LN(ABS("{score_column}") + 1e-300)))'

    con.execute(
        f"""
        CREATE TABLE prot_scores AS
        SELECT
            "{prot_col}",
            "{target_col}",
            {agg_expr} AS "{score_column}",
            random()     AS _crema_rand
        FROM unambig_psms
        GROUP BY "{prot_col}", "{target_col}"
        """
    )
    _compete_into(
        con, "prot_scores", "competed_prot", score_column, [prot_col], desc
    )
    prot_df = _fetch_qvalues(
        con, "competed_prot", score_column, target_col, desc
    )
    prot_df = prot_df.drop(columns=["_crema_rand"])
    mask_prot = prot_df[target_col].values.astype(bool)
    conf["proteins"] = prot_df[mask_prot].reset_index(drop=True)
    decoy_conf["proteins"] = prot_df[~mask_prot].reset_index(drop=True)
    LOGGER.info(
        "  - Found %i Proteins at q<=%g.",
        (conf["proteins"]["crema q-value"] <= eval_fdr).sum(),
        eval_fdr,
    )

    # Protein groups require a graph-based grouping algorithm that is not yet
    # implemented in the DuckDB backend.  The key is omitted so that the
    # writer and _prettify_tables simply skip it.

    con.close()
    return conf, decoy_conf
