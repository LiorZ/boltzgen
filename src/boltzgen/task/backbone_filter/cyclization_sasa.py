"""Cyclization filter: Cyclization site solvent exposure via rSASA."""

import warnings
from pathlib import Path

from Bio.PDB import MMCIFParser, ShrakeRupley

from .base import BaseFilter, FilterResult

# Max SASA per residue type (Tien et al. 2013, Gly-X-Gly tripeptide reference).
MAX_SASA = {
    "ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167,
    "GLN": 225, "GLU": 223, "GLY": 104, "HIS": 224, "ILE": 197,
    "LEU": 201, "LYS": 236, "MET": 224, "PHE": 240, "PRO": 159,
    "SER": 155, "THR": 172, "TRP": 285, "TYR": 263, "VAL": 174,
}


class CyclizationSolventExposureFilter(BaseFilter):
    """Verify that cyclization site terminal loops are solvent-exposed.

    Computes per-residue rSASA (relative SASA) for the terminal loop
    residues of the binder chain. The first residue (Cys) is excluded
    because it forms a disulfide in BoltzGen output.

    Args:
        binder_chain: Chain ID of the binder.
        n_terminal_loop: Number of loop residues to check at each terminus
            (excluding Cys-1 at the N-terminus).
        rsasa_threshold: Minimum rSASA for a residue to be considered exposed.
        min_exposed_count: Minimum number of exposed residues required at
            each terminus to pass (must be <= n_terminal_loop).
    """

    def __init__(
        self,
        binder_chain: str = "A",
        n_terminal_loop: int = 3,
        rsasa_threshold: float = 0.5,
        min_exposed_count: int = 2,
    ):
        self.binder_chain = binder_chain
        self.n_terminal_loop = n_terminal_loop
        self.rsasa_threshold = rsasa_threshold
        self.min_exposed_count = min_exposed_count

    @property
    def name(self) -> str:
        return "sasa"

    def run(self, cif_path: Path) -> FilterResult:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure("s", str(cif_path))
            model = structure[0]
            sr = ShrakeRupley()
            sr.compute(model, level="R")

        chain = model[self.binder_chain] if self.binder_chain in model else None
        if chain is None:
            return FilterResult(
                passed=False, metrics={}, details=f"Chain {self.binder_chain} not found"
            )

        residues = list(chain.get_residues())
        n_res = len(residues)
        if n_res < 2 * (self.n_terminal_loop + 1):
            return FilterResult(
                passed=False, metrics={}, details="Binder too short for SASA check"
            )

        # N-terminal loop: positions 2..2+n (skip Cys-1 at index 0)
        nterm_residues = residues[1 : 1 + self.n_terminal_loop]
        # C-terminal loop: last n residues
        cterm_residues = residues[-self.n_terminal_loop :]

        def _rsasa_values(res_list):
            values = []
            for res in res_list:
                sasa = res.sasa
                resname = res.get_resname()
                max_s = MAX_SASA.get(resname, 200)
                values.append(sasa / max_s)
            return values

        nterm_rsasa = _rsasa_values(nterm_residues)
        cterm_rsasa = _rsasa_values(cterm_residues)

        nterm_exposed = sum(1 for v in nterm_rsasa if v >= self.rsasa_threshold)
        cterm_exposed = sum(1 for v in cterm_rsasa if v >= self.rsasa_threshold)

        n_pass = nterm_exposed >= self.min_exposed_count
        c_pass = cterm_exposed >= self.min_exposed_count
        passed = n_pass and c_pass

        nterm_mean = sum(nterm_rsasa) / len(nterm_rsasa) if nterm_rsasa else 0.0
        cterm_mean = sum(cterm_rsasa) / len(cterm_rsasa) if cterm_rsasa else 0.0

        metrics = {
            "nterm_loop_mean_rsasa": round(nterm_mean, 3),
            "cterm_loop_mean_rsasa": round(cterm_mean, 3),
            "nterm_exposed_count": float(nterm_exposed),
            "cterm_exposed_count": float(cterm_exposed),
        }

        nterm_str = ", ".join(f"{v:.2f}" for v in nterm_rsasa)
        cterm_str = ", ".join(f"{v:.2f}" for v in cterm_rsasa)
        details = (
            f"N-loop rSASA: [{nterm_str}] ({nterm_exposed}/{len(nterm_rsasa)} exposed), "
            f"C-loop rSASA: [{cterm_str}] ({cterm_exposed}/{len(cterm_rsasa)} exposed)"
        )

        return FilterResult(passed=passed, metrics=metrics, details=details)
