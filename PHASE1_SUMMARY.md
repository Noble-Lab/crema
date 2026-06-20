# Phase 1 Implementation Summary

## What was done

Implemented **Phase 1** of `SCALING_PLAN.md` on branch `scale-phase-1-2`
(branched off `add-unit-tests`).  All 90 existing tests still pass.

---

## Changes made

### 1. Remove redundant copy from `PsmDataset.data` property
**File:** `crema/dataset.py`

- `data` property previously returned `self._data.copy()` on every access —
  wasting memory proportional to the full dataset on each call.
- Now returns `self._data` directly.
- Callers that genuinely need an independent copy must call `.data.copy()`
  themselves (none in the existing codebase do).

### 2. Change `copy_data` default from `True` to `False`
**Files:** `crema/dataset.py`, `crema/parsers/txt.py`, `crema/parsers/tide.py`,
`crema/parsers/comet.py`, `crema/parsers/msamanda.py`,
`crema/parsers/msfragger.py`, `crema/parsers/msgf.py`

- Every public entry-point (`PsmDataset.__init__`, `read_txt`, `read_tide`,
  `read_comet`, `read_msamanda`, `read_msfragger`, `read_msgf`) previously
  defaulted to making a full deep copy of any DataFrame passed in.
- Default is now `False`: no copy is made unless the caller explicitly passes
  `copy_data=True`.
- For the common case (reading from files), this has no effect because the
  DataFrame is freshly constructed from the CSV — there is nothing to copy.
- For callers passing an existing DataFrame, this saves one full copy.

### 3. Add optional pyarrow CSV engine
**Files:** `crema/parsers/txt.py` (`_parse_psms`), `crema/utils.py`
(`parse_psms_txt`)

- Both low-level CSV readers now try `engine="pyarrow"` first; if pyarrow is
  not installed they fall back to pandas' default C engine transparently.
- pyarrow typically parses 2–4× faster than the C engine on large files.
- No API change; pyarrow remains an optional dependency.

### 4. Minor test cleanup
**File:** `tests/unit_tests/test_confidence.py`

- Removed a now-stale comment that said "Use internal data to avoid the copy()
  on the .data property" — the copy no longer exists.

---

## How to push and create the PR

The commit is on local branch `scale-phase-1-2`.  Because SSH key access is
not available in this terminal session, push manually:

```bash
git push -u origin scale-phase-1-2
```

Then create the PR targeting `add-unit-tests`:

```bash
gh pr create \
  --base add-unit-tests \
  --title "Phase 1: eliminate redundant data copies and add optional pyarrow engine" \
  --body "$(cat <<'EOF'
## Summary

- Remove `.copy()` from `PsmDataset.data` property — previously allocated a
  full copy of the dataset on every property access.
- Change `copy_data` default to `False` in `PsmDataset.__init__` and all
  parser entry-points (`read_txt`, `read_tide`, `read_comet`, `read_msamanda`,
  `read_msfragger`, `read_msgf`) — eliminates a redundant deep copy when
  callers pass an existing DataFrame.
- Use pyarrow CSV engine (with C-engine fallback) in `_parse_psms()` and
  `parse_psms_txt()` — 2–4× faster parsing when pyarrow is installed.

All 90 unit tests pass.  This is Phase 1 of `SCALING_PLAN.md`.

## Test plan
- [ ] CI passes (90 unit tests + lint)
- [ ] Manually verify `read_tide` / `read_txt` still work on a real Tide file
- [ ] Install pyarrow and confirm pyarrow path is exercised (`LOGGER` shows no
  fallback warning)

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```
