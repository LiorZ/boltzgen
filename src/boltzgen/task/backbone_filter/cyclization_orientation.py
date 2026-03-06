"""Cyclization filter: Cyclization site orientation — CA->CB ray + distance to binding site."""

from pathlib import Path

import numpy as np

from .base import BaseFilter, FilterResult


def _parse_cif_atoms(cif_path: Path) -> dict[str, list[dict]]:
    """Parse CIF ATOM records into per-chain atom lists.

    Returns:
        Dict mapping chain_id to list of atom dicts with keys:
        atom_name, resnum, resname, x, y, z.
    """
    chains: dict[str, list[dict]] = {}
    with open(cif_path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 15 or parts[0] != "ATOM":
                continue
            atom_name = parts[3]
            chain = parts[6]
            resnum = int(parts[8])
            resname = parts[5]
            x, y, z = float(parts[10]), float(parts[11]), float(parts[12])
            chains.setdefault(chain, []).append(
                dict(atom_name=atom_name, resnum=resnum, resname=resname, x=x, y=y, z=z)
            )
    return chains


class CyclizationOrientationFilter(BaseFilter):
    """Reject scaffolds where cyclization site faces the target.

    Two sub-checks:
    1. **Ray test**: The CA->CB vector of the first and last binder residues
       is extended as a ray. If it passes within `ray_radius` of any target
       atom, the terminal is oriented toward the target -> fail.
    2. **Distance test**: Minimum CA-CA distance from terminal residues
       (first/last `n_terminal`) to target binding hotspot residues must
       exceed `min_distance`.

    Both sub-checks must pass.

    Args:
        binder_chain: Chain ID of the binder.
        target_chains: Chain ID(s) of the target.
        binding_residues: Target residue numbers defining the binding hotspot.
        n_terminal: Number of terminal residues for distance check.
        ray_radius: Maximum closest-approach distance (A) for the ray test.
        min_distance: Minimum CA-CA distance (A) to binding hotspot.
    """

    def __init__(
        self,
        binder_chain: str = "A",
        target_chains: list[str] | str | None = None,
        binding_residues: list[int] | None = None,
        n_terminal: int = 4,
        ray_radius: float = 8.0,
        min_distance: float = 10.0,
    ):
        self.binder_chain = binder_chain
        if target_chains is None:
            self.target_chains = ["B"]
        elif isinstance(target_chains, str):
            self.target_chains = [target_chains]
        else:
            self.target_chains = list(target_chains)
        self.binding_residues = set(binding_residues) if binding_residues else set()
        self.n_terminal = n_terminal
        self.ray_radius = ray_radius
        self.min_distance = min_distance

    @property
    def name(self) -> str:
        return "orientation"

    def _ray_min_approach(
        self, ca: np.ndarray, cb: np.ndarray, target_coords: np.ndarray
    ) -> float:
        """Minimum closest-approach distance from CA->CB ray to any target atom.

        The ray is: point = ca + t * (cb - ca), for t > 0 (forward only).
        """
        ray_dir = cb - ca
        ray_len = np.linalg.norm(ray_dir)
        if ray_len < 1e-6:
            return float("inf")
        ray_unit = ray_dir / ray_len

        to_target = target_coords - ca
        proj = np.dot(to_target, ray_unit)

        # Only forward direction (t > 0)
        forward_mask = proj > 0
        if not forward_mask.any():
            return float("inf")

        perp = to_target[forward_mask] - proj[forward_mask, None] * ray_unit
        dists = np.linalg.norm(perp, axis=1)
        return float(dists.min())

    def run(self, cif_path: Path) -> FilterResult:
        chains = _parse_cif_atoms(cif_path)

        if self.binder_chain not in chains:
            return FilterResult(
                passed=False, metrics={}, details=f"Binder chain {self.binder_chain} not found"
            )

        # Collect all target chain atoms
        target_atoms = []
        for tc in self.target_chains:
            if tc in chains:
                target_atoms.extend(chains[tc])
        if not target_atoms:
            return FilterResult(
                passed=False,
                metrics={},
                details=f"Target chain(s) {self.target_chains} not found",
            )

        binder_atoms = chains[self.binder_chain]

        # Build per-residue atom lookup for binder
        binder_by_res: dict[int, dict[str, np.ndarray]] = {}
        for a in binder_atoms:
            rn = a["resnum"]
            binder_by_res.setdefault(rn, {})[a["atom_name"]] = np.array(
                [a["x"], a["y"], a["z"]]
            )

        resnums = sorted(binder_by_res.keys())
        if len(resnums) < 2 * self.n_terminal:
            return FilterResult(
                passed=False, metrics={}, details="Binder too short for terminal check"
            )

        first_res = resnums[0]
        last_res = resnums[-1]

        # All target atom coords for ray test
        target_coords = np.array([[a["x"], a["y"], a["z"]] for a in target_atoms])

        # --- Ray test ---
        nterm_ray_dist = float("inf")
        cterm_ray_dist = float("inf")

        for res, label in [(first_res, "nterm"), (last_res, "cterm")]:
            atoms = binder_by_res[res]
            if "CA" in atoms and "CB" in atoms:
                d = self._ray_min_approach(atoms["CA"], atoms["CB"], target_coords)
                if label == "nterm":
                    nterm_ray_dist = d
                else:
                    cterm_ray_dist = d

        ray_hits = nterm_ray_dist < self.ray_radius or cterm_ray_dist < self.ray_radius

        # --- Distance test ---
        if self.binding_residues:
            hotspot_cas = []
            for a in target_atoms:
                if a["atom_name"] == "CA" and a["resnum"] in self.binding_residues:
                    hotspot_cas.append(np.array([a["x"], a["y"], a["z"]]))

            if not hotspot_cas:
                return FilterResult(
                    passed=False,
                    metrics={},
                    details=f"No binding hotspot CAs found for residues {self.binding_residues}",
                )

            hotspot_arr = np.array(hotspot_cas)
            term_resnums = resnums[: self.n_terminal] + resnums[-self.n_terminal :]

            min_term_to_binding = float("inf")
            for rn in term_resnums:
                if "CA" in binder_by_res[rn]:
                    ca = binder_by_res[rn]["CA"]
                    dists = np.linalg.norm(hotspot_arr - ca, axis=1)
                    d = float(dists.min())
                    if d < min_term_to_binding:
                        min_term_to_binding = d

            distance_pass = min_term_to_binding > self.min_distance
        else:
            # No binding residues specified — skip distance test
            min_term_to_binding = float("inf")
            distance_pass = True

        passed = (not ray_hits) and distance_pass

        metrics = {
            "nterm_ray_min_dist": round(nterm_ray_dist, 2),
            "cterm_ray_min_dist": round(cterm_ray_dist, 2),
            "ray_hits_target": float(ray_hits),
            "min_term_to_binding": round(min_term_to_binding, 2),
        }

        details = (
            f"Ray N-term: {nterm_ray_dist:.1f}A, C-term: {cterm_ray_dist:.1f}A "
            f"(radius={self.ray_radius}A, {'HIT' if ray_hits else 'clear'}); "
            f"Min dist to hotspot: {min_term_to_binding:.1f}A (thresh={self.min_distance}A)"
        )

        return FilterResult(passed=passed, metrics=metrics, details=details)
