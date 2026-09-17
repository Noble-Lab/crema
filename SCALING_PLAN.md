# Crema Scaling Plan: Supporting 20–30 Billion Target-Decoy Pairs

## Background

The current Crema codebase assumes all PSM data fits in RAM. At 50 bytes/PSM (minimal columns), 20–30 billion PSMs requires 1–1.5 TB of memory — far beyond what the current architecture can handle. This document describes a phased plan to scale Crema gracefully to datasets of that size.

---

## Current Architecture Bottlenecks

The following patterns in the existing code prevent scaling:

- **`PsmDataset._data`** is a single pandas DataFrame holding every PSM in memory.
- **`data` property** returns `self._data.copy()` — duplicating the entire dataset on every access.
- **Parsers** call `pd.concat([pd.read_csv(f) for f in files])` — all input files are loaded simultaneously into memory.
- **Competition phase** (`confidence.py::_compete`) performs an in-memory `groupby().apply()` over all PSMs.
- **`qvalues.tdc()`** calls `np.argsort()` on the full score array.

---

## Proposed Architecture

### Strategy: DuckDB as Out-of-Core Engine

[DuckDB](https://duckdb.org/) is the recommended engine for large-dataset operations. It:

- Handles out-of-core operations transparently (external sort, streaming window functions, parallel reads)
- Can read CSV, TSV, and Parquet files directly without loading them into Python
- Expresses the TDC algorithm naturally using SQL window functions
- Is a single Python dependency with no external server required
- Executes queries in parallel using all available cores by default

The plan has four phases, each independently useful.

---

## Phase 1 — Eliminate Unnecessary Copies

**Effort:** Low
**Impact:** 2–4× memory reduction; no API change
**Files:** `dataset.py`, `parsers/txt.py`, `confidence.py`

### Changes

1. **Remove the defensive copy from the `data` property.**
   `dataset.py` line 118 returns `self._data.copy()`. This duplicates the entire dataset on every property access. Change to return a view or the DataFrame directly. Callers that genuinely need a copy can call `.copy()` themselves.

2. **Default `copy_data=False` in all parsers.**
   All parsers currently default to `copy_data=True`, which creates a deep copy of the input DataFrame on construction. Change the default to `False`.

3. **Use memory-efficient dtypes throughout.**
   - `float32` instead of `float64` for score columns (sufficient precision for FDR ranking)
   - `bool` for the target/decoy column (currently stored as object in some paths)
   - `pd.Categorical` for peptide and protein string columns, which repeat heavily across PSMs

4. **Use the pyarrow CSV engine.**
   Replace `pd.read_csv(...)` with `pd.read_csv(..., engine='pyarrow')` for faster, lower-memory file I/O.

---

## Phase 2 — Streaming Parsers

**Effort:** Medium
**Impact:** Enables processing files larger than RAM during the parsing stage
**Files:** `parsers/*.py`, `dataset.py`, `writers/txt.py`

### Changes

1. **Add `chunk_size` parameter to all parsers.**
   When `chunk_size` is set, parsers yield chunks of PSMs rather than loading the whole file:

   ```python
   def read_txt(..., chunk_size=None):
       if chunk_size is None:
           # existing behavior (unchanged)
           ...
       else:
           for chunk in pd.read_csv(f, chunksize=chunk_size, usecols=cols):
               yield _make_psm_chunk(chunk)
   ```

2. **Switch XML parsers to iterative element parsing.**
   The pepXML and mzTab parsers load the entire XML document into memory. Switch to pyteomics' iterative parsing API, which processes one element at a time.

3. **Introduce `PsmDatasetStream`.**
   A new class (or a `lazy=True` mode on `PsmDataset`) that stores a DuckDB relation instead of a DataFrame. This becomes the internal representation for large datasets. The existing public API (`.data`, `.scores`, `.targets`, etc.) is preserved — properties query DuckDB instead of returning a DataFrame slice.

4. **Stream output from the writer.**
   `writers/txt.py` currently calls `df.to_csv()` on the full result DataFrame. Replace with streaming output from a DuckDB cursor, writing rows in batches.

---

## Phase 3 — DuckDB Backend for Competition and FDR

**Effort:** High
**Impact:** Enables end-to-end processing of datasets that do not fit in RAM
**Files:** `confidence.py`, `qvalues.py`, new `backends/duckdb.py`
**New dependency:** `duckdb`

This is the core change. The entire TDC pipeline maps cleanly to DuckDB SQL, which handles out-of-core execution automatically.

### Competition (best PSM per spectrum)

```sql
SELECT *
FROM (
  SELECT *,
    ROW_NUMBER() OVER (
      PARTITION BY file, scan
      ORDER BY score DESC
    ) AS rn
  FROM psms
)
WHERE rn = 1
```

DuckDB executes this with an external sort when the data exceeds available memory.

### TDC FDR + Q-values in one pass

```sql
WITH sorted AS (
  SELECT *,
    SUM(is_target::INT) OVER (
      ORDER BY score DESC ROWS UNBOUNDED PRECEDING
    ) AS cum_targets,
    SUM((NOT is_target)::INT) OVER (
      ORDER BY score DESC ROWS UNBOUNDED PRECEDING
    ) AS cum_decoys
  FROM competed
),
with_fdr AS (
  SELECT *,
    (cum_decoys + 1.0) / NULLIF(cum_targets, 0) AS fdr
  FROM sorted
)
SELECT *,
  -- backward monotone minimum = q-value
  MIN(fdr) OVER (
    ORDER BY score DESC
    ROWS BETWEEN CURRENT ROW AND UNBOUNDED FOLLOWING
  ) AS qvalue
FROM with_fdr
```

DuckDB's window functions handle the backward pass for q-values natively. This replaces `qvalues.tdc()` and `_fdr2qvalue()` for the large-data path. Tied scores require grouping before the window function, mirroring the existing logic in `qvalues.py:129–137`.

### Peptide-level FDR

Same SQL pattern: group competed PSMs by peptide sequence using target/decoy pairing (stored as a small DuckDB lookup table), then apply the window-function FDR query.

### Protein-level FDR

The current `_group_proteins()` function iterates in Python. Replace with SQL aggregation:

```sql
SELECT protein, MAX(peptide_qvalue) AS best_score, is_target
FROM peptide_protein_map
WHERE peptide_qvalue <= 0.01
GROUP BY protein, is_target
```

Then apply the same TDC window-function query to the protein-level rows.

### Mix-max confidence

The `mixmax` method involves bootstrap estimation of π₀ (`estimate_pi0()`), which is iterative. The bootstrap loop can remain in Python/NumPy, but the data feeding it should be sampled from DuckDB rather than materialized in full.

---

## Phase 4 — Parquet as Native Format

**Effort:** Low–Medium
**Impact:** 10× I/O speedup and ~10× storage reduction for repeat analyses
**Files:** New `parsers/parquet.py`, `writers/parquet.py`

Add Parquet as a first-class input/output format. DuckDB reads Parquet directly with zero-copy Arrow memory and predicate pushdown (only rows matching a filter are read from disk).

Recommended workflow for large datasets:

1. **First run:** Parse CSV/TSV input → write `psms.parquet` via `crema convert`
2. **Subsequent runs:** `crema assign-confidence psms.parquet` — reads directly, skips CSV parsing

This is particularly valuable for iterative analyses (e.g., varying the FDR threshold or score column) where re-parsing the raw files is the dominant cost.

---

## Summary of File Changes

| File | Change | Phase |
|---|---|---|
| `dataset.py` | Remove `.copy()` from property; add `PsmDatasetStream` backed by DuckDB | 1, 3 |
| `parsers/txt.py` | Add `chunk_size` param; use pyarrow engine; streaming yield mode | 1, 2 |
| `parsers/*.py` | Streaming mode; XML parsers switch to iterative element parsing | 2 |
| `qvalues.py` | Add `tdc_duckdb()` and `mixmax_duckdb()` variants; keep originals for small data | 3 |
| `confidence.py` | Route through DuckDB backend when dataset is large or `backend='duckdb'` | 3 |
| `writers/txt.py` | Stream output from DuckDB cursor in chunks | 2, 3 |
| `backends/duckdb.py` | New module: SQL queries for competition, FDR, and q-value calculation | 3 |
| `parsers/parquet.py` | New: read Parquet input | 4 |
| `writers/parquet.py` | New: write Parquet output | 4 |

---

## Backward Compatibility

- **Small datasets:** The existing pandas path remains the default. No API change for current users.
- **Large datasets:** Auto-switch to DuckDB when input exceeds a configurable threshold (e.g., `large_dataset_threshold=1_000_000` PSMs), or via explicit `backend='duckdb'` parameter on `assign_confidence()`.
- **Public API preserved:** `PsmDataset` properties (`.data`, `.scores`, `.targets`, etc.) continue to work the same way; they query DuckDB instead of slicing a DataFrame when in large-data mode.

---

## New Dependencies

| Package | Phase | Install |
|---|---|---|
| `pyarrow` | 1 | `pip install crema[fast]` |
| `duckdb` | 3 | `pip install crema[large]` |

Phases 1 and 2 require no new mandatory dependencies. Phase 3 adds `duckdb` as an optional dependency. Both packages ship as pure Python wheels with no C compiler required.
