"""Fixed-score corrections and reference-conditioned leakage diagnostics."""
from .core import analyze_bulk_scores, load_paper_reference, load_paper_signatures

__version__ = "0.1.0"
__all__ = ["analyze_bulk_scores", "load_paper_reference", "load_paper_signatures"]
