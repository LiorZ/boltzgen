"""Cyclization filter: Terminal Loop Verification via DSSP secondary structure assignment."""

import warnings
from pathlib import Path

import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx

from .base import BaseFilter, FilterResult


class CyclizationTerminalLoopFilter(BaseFilter):
    """Verify that binder terminal residues are in loop/coil conformation.

    Uses biotite's annotate_sse() (DSSP-based) to assign secondary structure.
    Checks that the first and last N residues of the binder chain are
    predominantly coil ('c'), not helix ('a') or strand ('b').

    Args:
        binder_chain: Chain ID of the binder in the CIF file.
        n_terminal: Number of residues to check at each terminus.
        min_coil_count: Minimum number of coil residues required at each
            terminus to pass (must be <= n_terminal).
    """

    def __init__(
        self,
        binder_chain: str = "A",
        n_terminal: int = 4,
        min_coil_count: int = 3,
    ):
        self.binder_chain = binder_chain
        self.n_terminal = n_terminal
        self.min_coil_count = min_coil_count

    @property
    def name(self) -> str:
        return "terminal_ss"

    def run(self, cif_path: Path) -> FilterResult:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cif_file = pdbx.CIFFile.read(str(cif_path))
            structure = pdbx.get_structure(cif_file, model=1)

        binder = structure[structure.chain_id == self.binder_chain]
        if len(binder) == 0:
            return FilterResult(
                passed=False, metrics={}, details=f"Chain {self.binder_chain} not found"
            )

        sse = struc.annotate_sse(binder)  # one entry per residue
        n = self.n_terminal

        n_term_sse = sse[:n]
        c_term_sse = sse[-n:]

        n_term_coil = int(sum(1 for s in n_term_sse if s == "c"))
        n_term_helix = int(sum(1 for s in n_term_sse if s == "a"))
        n_term_strand = int(sum(1 for s in n_term_sse if s == "b"))

        c_term_coil = int(sum(1 for s in c_term_sse if s == "c"))
        c_term_helix = int(sum(1 for s in c_term_sse if s == "a"))
        c_term_strand = int(sum(1 for s in c_term_sse if s == "b"))

        n_pass = n_term_coil >= self.min_coil_count
        c_pass = c_term_coil >= self.min_coil_count
        passed = n_pass and c_pass

        metrics = {
            "n_term_coil_count": n_term_coil,
            "n_term_helix_count": n_term_helix,
            "n_term_strand_count": n_term_strand,
            "c_term_coil_count": c_term_coil,
            "c_term_helix_count": c_term_helix,
            "c_term_strand_count": c_term_strand,
        }

        details = (
            f"N-term SSE: {''.join(n_term_sse)} (coil={n_term_coil}/{n}), "
            f"C-term SSE: {''.join(c_term_sse)} (coil={c_term_coil}/{n})"
        )

        return FilterResult(passed=passed, metrics=metrics, details=details)
