# Agentic Roadmap for Interactions-search

This roadmap is designed for **small, independent agent sessions** (Claude Code / Copilot / similar), each one with clear scope and acceptance criteria.

## Current decisions

- Product direction: **CLI and importable library are both first-class**
- Python baseline: **3.11+**
- Language: **English**
- Config approach: **YAML remains default**, CLI flags can override values

---

## How to run this roadmap

For each step:
1. Open a new agent session.
2. Paste the "Agent prompt" block for that step.
3. Ask the agent to stop after validation and summary.
4. Review diff before moving to next step.

Keep each step small. If a step feels too large, split it.

---

## Step 1 — Baseline snapshot (no refactor)

**Goal**
Capture current behavior to avoid regressions during refactors.

**Files to create**
- `tests/smoke/test_cli_invocation.py`
- `tests/data/` with 1 minimal receptor + ligand fixture pair
- `docs/BASELINE.md`

**Acceptance criteria**
- Running current script on fixture data succeeds.
- Output files are generated as expected.
- Baseline command and expected artifacts documented.

**Agent prompt**
```text
Create a minimal smoke test and baseline docs for the current Interactions_search.py behavior.
Do not refactor production code yet.
Add tiny test fixtures and verify CLI execution on them.
Document command, outputs, and known limitations in docs/BASELINE.md.
```

---

## Step 2 — Package skeleton (still compatible)

**Goal**
Create importable package structure while preserving current script behavior.

**Files to create**
- `src/interactions_search/__init__.py`
- `src/interactions_search/cli.py`
- `src/interactions_search/config.py`
- `src/interactions_search/types.py`

**Acceptance criteria**
- `python Interactions_search.py ...` still works.
- `python -m interactions_search.cli --help` works.
- No logic rewrite; mostly moving wrappers and interfaces.

**Agent prompt**
```text
Introduce src/interactions_search package with cli/config/types modules.
Keep old Interactions_search.py path working.
Do minimal adapters only; no algorithmic changes.
Wire argparse into a callable main() in the new package.
```

---

## Step 3 — Typed config + YAML override model

**Goal**
Replace loose dict config with a typed config object and deterministic overrides.

**Files to create/update**
- `src/interactions_search/config.py`
- `tests/config/test_config_loading.py`
- `docs/CONFIG.md`

**Acceptance criteria**
- YAML loads into typed model.
- CLI flags override YAML keys predictably.
- Invalid config fails with clear errors.

**Agent prompt**
```text
Implement typed configuration loading from Interacciones_variables.yml, with optional CLI overrides.
Add unit tests for precedence rules and invalid values.
Document final precedence in docs/CONFIG.md.
```

---

## Step 4 — Core engine extraction

**Goal**
Move domain logic from monolith script to focused modules.

**Target modules**
- `src/interactions_search/io_pdb.py`
- `src/interactions_search/geometry.py`
- `src/interactions_search/interactions.py`
- `src/interactions_search/reporting.py`

**Acceptance criteria**
- Same outputs for baseline fixture.
- Functions are importable and testable independently.
- No silent behavior changes.

**Agent prompt**
```text
Refactor Interactions_search.py into focused modules under src/interactions_search.
Preserve behavior verified in baseline tests.
Prefer moving code first, then tiny cleanups.
Add or update tests only where needed to prove parity.
```

---

## Step 5 — Quality gates + CI

**Goal**
Automate checks so every future change is safer.

**Files to create**
- `.github/workflows/ci.yml`
- `tests/` expansion for core functions

**Acceptance criteria**
- CI runs lint + tests on Python 3.11.
- Failing checks block regressions.

**Agent prompt**
```text
Add GitHub Actions CI for pytest and ruff.
Keep matrix minimal (start with Python 3.11).
Make sure project installs from pyproject.toml in CI.
```

---

## Step 6 — Library API + release readiness

**Goal**
Expose stable import API and prepare versioned releases.

**Files to create/update**
- `src/interactions_search/api.py`
- `src/interactions_search/__init__.py`
- `docs/API.md`
- `CHANGELOG.md`

**Acceptance criteria**
- Public API documented and stable.
- CLI and library workflows both covered in docs.

**Agent prompt**
```text
Define a minimal stable library API in interactions_search.api and export it in __init__.
Document usage examples for both CLI and Python import usage.
Create initial changelog with semantic versioning notes.
```

---

## Guardrails for every agent step

- Keep changes scoped to one step.
- Run tests/lint after edits.
- Do not mix refactor + feature additions in one PR.
- Preserve scientific criteria unless explicitly changed and validated.
- Update docs together with code.

