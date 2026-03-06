"""BackboneFilter Task — pipeline step that filters generated backbones before inverse folding."""

import shutil
from pathlib import Path
from typing import List, Optional

from omegaconf import OmegaConf

from boltzgen.task.task import Task
from .base import BaseFilter
from .cyclization_orientation import CyclizationOrientationFilter
from .pipeline import FilterPipeline
from .cyclization_sasa import CyclizationSolventExposureFilter
from .cyclization_terminal_ss import CyclizationTerminalLoopFilter


class BackboneFilter(Task):
    """Filter generated backbone CIFs between design and inverse folding steps.

    Runs a configurable set of structural filters on CIF files in `input_dir`.
    Failing structures are moved to a subfolder so downstream steps skip them.

    All parameters are fully configurable via YAML / CLI — no hardcoded
    residue indices, chain IDs, or thresholds.

    Parameters
    ----------
    input_dir : str
        Directory containing CIF files from the design step.
    binder_chain : str
        Chain ID of the designed binder in the CIF files.
    target_chains : list[str]
        Chain ID(s) of the target protein.
    terminal_ss_enabled : bool
        Enable the terminal secondary structure filter.
    ss_n_terminal : int
        Number of residues to check at each terminus for SS.
    ss_min_coil : int
        Minimum coil residues required at each terminus.
    orientation_enabled : bool
        Enable the orientation filter (ray + distance).
    binding_residues : list[int] or None
        Target residue numbers defining the binding hotspot.
    ray_radius : float
        CA->CB ray proximity threshold (Angstroms).
    min_distance : float
        Minimum terminal-to-hotspot CA distance (Angstroms).
    orient_n_terminal : int
        Number of terminal residues for the distance check.
    sasa_enabled : bool
        Enable the solvent exposure filter.
    rsasa_threshold : float
        Minimum rSASA for a residue to count as exposed.
    sasa_n_loop : int
        Number of loop residues to check at each terminus.
    sasa_min_exposed : int
        Minimum exposed residues required at each terminus.
    move_failed : bool
        If True, move failing CIFs to a subfolder instead of deleting.
    glob_pattern : str
        Glob pattern to match CIF files.
    """

    def __init__(
        self,
        input_dir: str,
        # Chain identification
        binder_chain: str = "A",
        target_chains: Optional[List[str]] = None,
        # Terminal SS filter
        terminal_ss_enabled: bool = True,
        ss_n_terminal: int = 4,
        ss_min_coil: int = 3,
        # Orientation filter
        orientation_enabled: bool = True,
        binding_residues: Optional[List[int]] = None,
        ray_radius: float = 8.0,
        min_distance: float = 10.0,
        orient_n_terminal: int = 4,
        # SASA filter
        sasa_enabled: bool = True,
        rsasa_threshold: float = 0.5,
        sasa_n_loop: int = 3,
        sasa_min_exposed: int = 2,
        # Behavior
        move_failed: bool = True,
        glob_pattern: str = "*.cif",
    ):
        self.input_dir = Path(input_dir)
        self.binder_chain = binder_chain
        self.target_chains = target_chains or ["B"]
        self.terminal_ss_enabled = terminal_ss_enabled
        self.ss_n_terminal = ss_n_terminal
        self.ss_min_coil = ss_min_coil
        self.orientation_enabled = orientation_enabled
        self.binding_residues = binding_residues
        self.ray_radius = ray_radius
        self.min_distance = min_distance
        self.orient_n_terminal = orient_n_terminal
        self.sasa_enabled = sasa_enabled
        self.rsasa_threshold = rsasa_threshold
        self.sasa_n_loop = sasa_n_loop
        self.sasa_min_exposed = sasa_min_exposed
        self.move_failed = move_failed
        self.glob_pattern = glob_pattern

    def _build_filters(self) -> list[BaseFilter]:
        """Build the list of enabled filters from config."""
        filters: list[BaseFilter] = []

        if self.terminal_ss_enabled:
            filters.append(
                CyclizationTerminalLoopFilter(
                    binder_chain=self.binder_chain,
                    n_terminal=self.ss_n_terminal,
                    min_coil_count=self.ss_min_coil,
                )
            )

        if self.orientation_enabled:
            filters.append(
                CyclizationOrientationFilter(
                    binder_chain=self.binder_chain,
                    target_chains=self.target_chains,
                    binding_residues=self.binding_residues,
                    n_terminal=self.orient_n_terminal,
                    ray_radius=self.ray_radius,
                    min_distance=self.min_distance,
                )
            )

        if self.sasa_enabled:
            filters.append(
                CyclizationSolventExposureFilter(
                    binder_chain=self.binder_chain,
                    n_terminal_loop=self.sasa_n_loop,
                    rsasa_threshold=self.rsasa_threshold,
                    min_exposed_count=self.sasa_min_exposed,
                )
            )

        return filters

    def run(self, config: OmegaConf = None) -> None:
        """Run the backbone filter pipeline.

        1. Build filters from config.
        2. Run on all CIFs in input_dir.
        3. Move/remove failing CIFs so inverse folding skips them.
        """
        filters = self._build_filters()
        if not filters:
            print("No backbone filters enabled — skipping.")
            return

        print(f"\nBackbone filtering: {len(filters)} filter(s) enabled")
        for f in filters:
            print(f"  - {f.name}")

        pipeline = FilterPipeline(filters)
        output_csv = self.input_dir / "backbone_filter_results.csv"

        df = pipeline.run_directory(
            input_dir=self.input_dir,
            output_csv=output_csv,
            glob_pattern=self.glob_pattern,
        )

        if df.empty:
            print("No structures found to filter.")
            return

        # Move failing CIFs
        failed_files = df[~df["all_pass"]]["file"].tolist()
        passed_files = df[df["all_pass"]]["file"].tolist()

        if failed_files and self.move_failed:
            reject_dir = self.input_dir / "backbone_filtered_out"
            reject_dir.mkdir(exist_ok=True)

            for fname in failed_files:
                src = self.input_dir / fname
                if src.exists():
                    shutil.move(str(src), str(reject_dir / fname))
                # Also move associated metadata files
                npz = src.with_suffix(".npz")
                if npz.exists():
                    shutil.move(str(npz), str(reject_dir / npz.name))
                native = self.input_dir / f"{src.stem}_native.cif"
                if native.exists():
                    shutil.move(str(native), str(reject_dir / native.name))

            print(
                f"\nMoved {len(failed_files)} failing backbone(s) to {reject_dir}"
            )

        print(
            f"Backbone filter complete: {len(passed_files)} passed, "
            f"{len(failed_files)} failed out of {len(df)} total"
        )
