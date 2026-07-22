import sys
from pathlib import Path

# Interactions_search.py lives at the repo root, two levels above src/interactions_search/
_repo_root = Path(__file__).parent.parent.parent
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))

from Interactions_search import main, analyze_pair, carga_variables  # noqa: E402

__all__ = ["main", "analyze_pair", "carga_variables"]
