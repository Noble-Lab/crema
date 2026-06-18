# Crema Unit Testing Plan

## Current State

The test suite has four unit test files and two system test files. Here is what is already covered and what is not.

### What Already Exists

| File | What It Tests | Quality |
|---|---|---|
| `test_qvalues.py` | `tdc()` and `mixmax()` numerical correctness, basic input validation | Good |
| `test_dataset.py` | `PsmDataset` construction, properties, `__getitem__`, `find_best_score()` | Partial |
| `test_parsers.py` | Each parser format (Tide, MSGF+, MSAmanda, MSFragger, Comet, mzTab, pepXML) | Good |
| `test_confidence.py` | `assign_confidence()` returns correct class type — **no assertions on values** | Thin |
| `test_cli.py` (system) | CLI invocation, file existence in output | Basic |
| `test_crux_tdc.py` (system) | Crux comparison | Unknown coverage |

### Key Gaps

1. **`test_confidence.py`** has `# TODO: assertions` comments — no numerical correctness tests at all for the FDR pipeline output.
2. **Multi-level FDR** (peptide, protein, protein-group) is entirely untested.
3. **Writers** (`writers/txt.py`) have no tests.
4. **Edge cases** (empty data, single PSM, all targets win, all decoys win, tied scores) are sparse.
5. **Error handling** (missing columns, malformed files, invalid arguments) is minimally tested.
6. **Multiple file inputs** and concatenation behavior are untested.

---

## Plan

### Infrastructure First

Before adding tests, add two pytest plugins to `pyproject.toml` (or `setup.cfg`):

```toml
[project.optional-dependencies]
dev = [
    "pytest",
    "pytest-cov",     # coverage reporting
    "hypothesis",     # property-based testing for numerical algorithms
]
```

Run tests with coverage via:
```
pytest --cov=crema --cov-report=term-missing tests/
```

---

### 1. `test_qvalues.py` — Fill Edge Cases

**What to add:**

#### 1a. TDC edge cases (parametrize these together)

| Scenario | Input | Expected behavior |
|---|---|---|
| All targets win | scores=[3,2,1], targets=[T,T,T] | All q-values = 1/N |
| All decoys win | scores=[3,2,1], targets=[D,D,D] | q-value = 1.0 for all |
| Single target PSM | scores=[1], targets=[T] | q-value = 1.0 |
| Single decoy PSM | scores=[1], targets=[D] | q-value = 1.0 |
| Two PSMs, target beats decoy | scores=[2,1], targets=[T,D] | q-value = 0.5, 0.5 |
| Two PSMs, decoy beats target | scores=[2,1], targets=[D,T] | q-value = 1.0, 1.0 |
| All tied scores | scores=[5,5,5,5], targets=[T,T,D,D] | Verify grouping logic |
| Large tied block at top | Prepend 100 tied targets to existing fixture | First block q-value ≤ existing min |

Use `@pytest.mark.parametrize` to cover all of these systematically.

#### 1b. Numerical stability

- Scores of `float('inf')` and `-float('inf')` — should not crash
- Scores of `float('nan')` — document and test expected behavior (likely undefined; at minimum should not silently give wrong answers)
- Very large scores (e.g., `1e308`) — should not overflow

#### 1c. FDR monotonicity invariant

Add a property-based test using `hypothesis`:

```python
from hypothesis import given, strategies as st

@given(st.lists(st.floats(allow_nan=False, allow_infinity=False), min_size=2),
       st.lists(st.booleans(), min_size=2))
def test_tdc_qvalues_monotone(scores, targets):
    """Q-values must be non-decreasing as score decreases."""
    # Align lengths, ensure at least one target
    ...
    qvals = tdc(np.array(scores[:n]), np.array(targets[:n]), desc=True)
    sorted_qvals = qvals[np.argsort(-np.array(scores[:n]))]
    assert all(sorted_qvals[i] <= sorted_qvals[i+1] for i in range(len(sorted_qvals)-1))
```

This catches regressions in the monotone minimum logic that hand-crafted examples might miss.

#### 1d. `_fdr2qvalue` directly

The Numba-compiled `_fdr2qvalue` is called internally but never tested directly. Add a test that imports it and verifies it produces a monotone non-increasing sequence given a non-monotone FDR array.

---

### 2. `test_dataset.py` — Validation and Error Handling

**What to add:**

#### 2a. Construction errors

| Scenario | Expected |
|---|---|
| `target_column` not in DataFrame | `ValueError` or `KeyError` |
| `score_columns` contains a column not in DataFrame | `ValueError` |
| `spectrum_columns` not in DataFrame | `ValueError` |
| Empty DataFrame | Raises or returns empty `PsmDataset` (document which) |
| Target column with values other than True/False/0/1 | `ValueError` |

#### 2b. Counts

```python
def test_target_decoy_counts(simple_psms):
    assert simple_psms._num_targets == 6
    assert simple_psms._num_decoys == 4
```

These are used in FDR calculation and should be tested explicitly.

#### 2c. `peptide_pairing` property

Currently `PsmDataset` holds `_peptide_pairing` but it is only tested indirectly through `read_tide()`. Add a direct test that sets it via the constructor and verifies it is returned correctly.

#### 2d. `find_best_score()` edge cases

- When all score columns perform equally (tie-breaking behavior)
- When `eval_fdr=0.0` (zero PSMs pass)
- When `eval_fdr=1.0` (all PSMs pass)

---

### 3. `test_parsers.py` — Error Handling and Multi-File Inputs

**What to add:**

#### 3a. Multiple file inputs

```python
def test_read_tide_multiple_files(target_tide_txt, decoy_tide_txt):
    """Concatenating two files should give combined PSM count."""
    psms = crema.read_tide([target_tide_txt, decoy_tide_txt])
    assert psms.data.shape[0] == 20  # 10 + 10
    assert psms.targets.sum() == 10
    assert (~psms.targets).sum() == 10
```

Test this for at least `read_tide` and `read_txt`. This exercises the `pd.concat` path that currently has no test.

#### 3b. Missing required columns

```python
def test_read_txt_missing_score_column(basic_tide_csv, tmp_path):
    with pytest.raises((KeyError, ValueError)):
        crema.read_txt(basic_tide_csv, ..., score_columns=["nonexistent"])
```

#### 3c. `copy_data=False`

Verify that `copy_data=False` returns a `PsmDataset` with the same data (not a copy). This is important for the Phase 1 memory optimization in the scaling plan.

```python
def test_read_txt_no_copy(basic_tide_df, tmp_path):
    ...
    psms = crema.read_txt(file, ..., copy_data=False)
    assert isinstance(psms.data, pd.DataFrame)
    assert psms.data.shape == (10, 6)
```

#### 3d. `protein_delim` parsing

Verify that proteins with delimiters (e.g., `"p1,p2"`) are correctly split or preserved depending on context, and that `protein_delim` is stored correctly.

#### 3e. Tide-specific: only targets or only decoys file

```python
def test_read_tide_target_only(target_tide_txt):
    """Single-file Tide input with only targets should work."""
    psms = crema.read_tide(target_tide_txt)
    assert psms.targets.all()
```

---

### 4. `test_confidence.py` — Numerical Correctness (Highest Priority Gap)

This is the most important section to add. The current tests only check that the return type is correct.

**What to add:**

#### 4a. Hand-computed TDC correctness at PSM level

Using `basic_tide_df` from `conftest.py` (10 PSMs, known scores), compute expected q-values by hand and assert them:

```python
def test_tdc_psm_qvalues(simple_psms):
    conf = simple_psms.assign_confidence(score_column="x", method="tdc")
    psm_results = conf.confidence_estimates["psms"]

    # Filter to targets
    targets = psm_results[psm_results["crema q-value"] <= 0.5]
    assert len(targets) == EXPECTED_COUNT

    # Check specific q-values
    np.testing.assert_allclose(
        psm_results.sort_values("x", ascending=False)["crema q-value"].values,
        EXPECTED_QVALUES,
        atol=1e-6
    )
```

Compute `EXPECTED_QVALUES` by hand from the fixture data (the fixture has known scores and target/decoy labels, so the TDC result is deterministic).

#### 4b. `pep_fdr_type` options

```python
@pytest.mark.parametrize("pep_fdr_type", ["psm-only", "peptide-only", "psm-peptide"])
def test_pep_fdr_type(simple_psms_with_pairing, pep_fdr_type):
    conf = simple_psms_with_pairing.assign_confidence(
        score_column="x", method="tdc", pep_fdr_type=pep_fdr_type
    )
    assert "peptides" in conf.confidence_estimates
    # More specific assertions per type
```

#### 4c. Peptide-level FDR correctness

Peptide-level results should:
- Have at most one row per unique target peptide
- Have q-values that are non-decreasing as score decreases
- Have `_num_peptides` that is ≤ `_num_targets`

#### 4d. Protein-level FDR correctness

```python
def test_protein_level_fdr(simple_psms):
    conf = simple_psms.assign_confidence(score_column="x", method="tdc")
    protein_results = conf.confidence_estimates["proteins"]

    # Each protein should appear at most once
    assert protein_results["protein id"].nunique() == len(protein_results)

    # Q-values should be in [0, 1]
    assert (protein_results["crema q-value"] >= 0).all()
    assert (protein_results["crema q-value"] <= 1).all()
```

#### 4e. `protein_score` aggregation methods

```python
@pytest.mark.parametrize("protein_score", ["max", "min", "sum", "prod"])
def test_protein_score_aggregation(simple_psms, protein_score):
    conf = simple_psms.assign_confidence(
        score_column="x", method="tdc", protein_score=protein_score
    )
    assert "proteins" in conf.confidence_estimates
```

#### 4f. Threshold filtering

```python
@pytest.mark.parametrize("fdr", [0.01, 0.05, 0.1, 0.5])
def test_fdr_threshold(simple_psms, fdr):
    conf = simple_psms.assign_confidence(score_column="x", method="tdc", eval_fdr=fdr)
    psm_results = conf.confidence_estimates["psms"]
    passing = psm_results[psm_results["crema q-value"] <= fdr]
    # All passing PSMs should actually be targets
    assert passing["target"].all()
```

#### 4g. Monotonicity invariant across all levels

For any confidence result, q-values must be non-decreasing as score decreases. Add a helper and call it for PSM, peptide, and protein levels.

---

### 5. New File: `test_writers.py`

The output writer (`writers/txt.py`) has no tests at all.

**What to add:**

#### 5a. Output file content

```python
def test_to_txt_creates_files(simple_psms, tmp_path):
    conf = simple_psms.assign_confidence(score_column="x")
    conf.to_txt(output_dir=tmp_path, file_root="test")

    assert (tmp_path / "test.crema.psms.txt").exists()
    assert (tmp_path / "test.crema.peptides.txt").exists()
    assert (tmp_path / "test.crema.proteins.txt").exists()
```

#### 5b. Output file column names

```python
def test_to_txt_columns(simple_psms, tmp_path):
    conf = simple_psms.assign_confidence(score_column="x")
    conf.to_txt(output_dir=tmp_path, file_root="test")

    result = pd.read_csv(tmp_path / "test.crema.psms.txt", sep="\t")
    assert "crema q-value" in result.columns
    assert "x" in result.columns
```

#### 5c. Round-trip: write then read

Confirm that q-values written to a file, read back, and compared to the in-memory values are identical (no floating-point corruption from string formatting).

#### 5d. FDR threshold filtering in output

When a `fdr_threshold` is passed to `to_txt()`, confirm only PSMs with q-value ≤ threshold appear in the output.

---

### 6. New File: `test_edge_cases.py`

Collect edge cases that cut across multiple modules.

#### 6a. Empty dataset

```python
def test_empty_dataset():
    """An empty PsmDataset should not crash FDR calculation."""
    empty_df = pd.DataFrame(columns=["scan", "score", "target", "peptide", "protein"])
    psms = PsmDataset(psms=empty_df, ...)
    # assign_confidence should either raise a clear error or return empty results
```

#### 6b. Single spectrum with one target and one decoy

The most minimal meaningful case. The target wins: one PSM passes with q-value = 1.0. Verify this explicitly.

#### 6c. Multiple PSMs per spectrum — only best survives competition

```python
def test_competition_keeps_best(tmp_path):
    """When multiple PSMs compete for one spectrum, only the best score survives."""
    # Build DataFrame with 3 PSMs for scan=1: scores 0.9 (decoy), 0.8 (target), 0.1 (target)
    # After competition, the decoy at 0.9 should win scan=1
    ...
    conf = psms.assign_confidence(score_column="x", method="tdc")
    psm_results = conf.confidence_estimates["psms"]
    scan1 = psm_results[psm_results["scan"] == 1]
    assert len(scan1) == 1
    assert not scan1["target"].iloc[0]  # decoy won
```

#### 6d. Protein grouping with shared peptides

A peptide that maps to multiple proteins should be excluded from protein-level FDR per the current implementation. Verify this:

```python
def test_shared_peptides_excluded_from_proteins():
    # PSM for peptide "APPLE" maps to protein "p1,p2" (protein_delim=",")
    # This peptide should not contribute to protein-level FDR
    ...
```

---

## Summary of New Test Count

| File | New Tests | Priority |
|---|---|---|
| `test_qvalues.py` | ~10 (edge cases + 1 property-based) | Medium |
| `test_dataset.py` | ~8 (validation + counts) | Medium |
| `test_parsers.py` | ~6 (multi-file, error handling) | Medium |
| `test_confidence.py` | ~12 (numerical correctness, all levels) | **High** |
| `test_writers.py` (new) | ~5 | Medium |
| `test_edge_cases.py` (new) | ~6 | Medium |
| **Total** | **~47** | |

## Recommended Implementation Order

1. **`test_confidence.py` numerical correctness** — This is the most important gap. The FDR pipeline output is the core deliverable of Crema and is currently tested only for return type.
2. **`test_writers.py`** — A natural companion to confidence tests; confirms the full pipeline end-to-end.
3. **`test_edge_cases.py`** — Catches regressions in corner cases that are easy to break during the scaling refactor.
4. **`test_qvalues.py` additions** — The existing tests are good; these additions harden the numerical core.
5. **`test_dataset.py` and `test_parsers.py` additions** — Important for robustness but lower risk area.
