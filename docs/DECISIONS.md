# Architecture Decisions Log

## 2026-07-20

### Decision 1: Product shape
- We support **both**:
  - CLI workflows for batch analysis
  - Python library usage for integration in other tools/pipelines

### Decision 2: Python compatibility
- Minimum supported version: **Python 3.11**

### Decision 3: Configuration strategy
- Keep YAML as the canonical source for defaults.
- Add CLI flags to override selected YAML values at runtime.
- Define and document deterministic precedence rules.

### Decision 4: Language
- Code comments and documentation will move toward **English-only**.
