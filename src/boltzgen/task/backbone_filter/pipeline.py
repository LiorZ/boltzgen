"""FilterPipeline — run all filters on a directory of CIFs and aggregate results."""

import csv
import time
from pathlib import Path

import pandas as pd

from .base import BaseFilter, FilterResult


class FilterPipeline:
    """Run a list of filters on all CIF files in a directory.

    Args:
        filters: List of BaseFilter instances to apply.
    """

    def __init__(self, filters: list[BaseFilter]):
        self.filters = filters

    def run_single(self, cif_path: Path) -> dict:
        """Run all filters on one CIF and return a flat dict of results."""
        row: dict = {"file": cif_path.name}

        all_passed = True
        for filt in self.filters:
            result = filt.run(cif_path)
            prefix = filt.name
            row[f"{prefix}_pass"] = result.passed
            row[f"{prefix}_details"] = result.details
            for k, v in result.metrics.items():
                row[f"{prefix}_{k}"] = v
            if not result.passed:
                all_passed = False

        row["all_pass"] = all_passed
        return row

    def run_directory(
        self,
        input_dir: Path,
        output_csv: Path | None = None,
        glob_pattern: str = "*.cif",
    ) -> pd.DataFrame:
        """Run all filters on every CIF in a directory.

        Writes results incrementally to CSV so progress is not lost on
        interruption. Prints progress every 100 structures.

        Args:
            input_dir: Directory containing CIF files.
            output_csv: If provided, write results to this CSV path.
            glob_pattern: Glob pattern to match CIF files.

        Returns:
            DataFrame with one row per CIF and columns for each filter's
            pass/fail, metrics, and details.
        """
        cif_files = sorted(input_dir.glob(glob_pattern))
        if not cif_files:
            print(f"No CIF files found in {input_dir} matching '{glob_pattern}'")
            return pd.DataFrame()

        total = len(cif_files)
        rows = []
        csv_writer = None
        csv_file = None
        header_written = False
        t0 = time.time()

        try:
            for i, cif_path in enumerate(cif_files):
                row = self.run_single(cif_path)
                rows.append(row)

                # Incremental CSV write
                if output_csv is not None:
                    if not header_written:
                        output_csv.parent.mkdir(parents=True, exist_ok=True)
                        csv_file = open(output_csv, "w", newline="")
                        csv_writer = csv.DictWriter(csv_file, fieldnames=row.keys())
                        csv_writer.writeheader()
                        header_written = True
                    csv_writer.writerow(row)
                    csv_file.flush()

                # Progress
                if (i + 1) % 100 == 0 or (i + 1) == total:
                    elapsed = time.time() - t0
                    rate = (i + 1) / elapsed
                    eta = (total - i - 1) / rate if rate > 0 else 0
                    n_pass = sum(1 for r in rows if r["all_pass"])
                    print(
                        f"  [{i+1:>6d}/{total}] "
                        f"{elapsed:.0f}s elapsed, {eta:.0f}s remaining, "
                        f"{rate:.1f} struct/s | "
                        f"{n_pass} passing all filters",
                        flush=True,
                    )
        finally:
            if csv_file is not None:
                csv_file.close()

        if output_csv is not None:
            print(f"\nResults written to {output_csv}")

        df = pd.DataFrame(rows)
        self._print_summary(df)
        return df

    def _print_summary(self, df: pd.DataFrame) -> None:
        """Print filter pass rates and overall statistics."""
        n = len(df)
        print(f"\n{'=' * 60}")
        print(f"Backbone Filter Summary — {n} structures")
        print(f"{'=' * 60}")

        for filt in self.filters:
            col = f"{filt.name}_pass"
            if col in df.columns:
                n_pass = df[col].sum()
                print(f"  {filt.name:20s}: {int(n_pass):4d} / {n} pass ({100 * n_pass / n:.1f}%)")

        if "all_pass" in df.columns:
            n_all = df["all_pass"].sum()
            print(f"  {'ALL FILTERS':20s}: {int(n_all):4d} / {n} pass ({100 * n_all / n:.1f}%)")

        print(f"{'=' * 60}")
