# Complex RMSD Filtering Based on Diffusion Samples

This document explains how to use the complex RMSD filtering feature, which enables filtering of refolded complexes based on the number of diffusion samples that pass a specified RMSD threshold.

## Overview

The complex RMSD filtering feature verifies that the binding conformation is recapitulated across multiple diffusion samples. This provides more robust validation than using only the best sample.

### How It Works

The RMSD calculation differs from standard RMSD in two key ways:

1. **Alignment**: The non-designed chains (receptor/target protein) are aligned first using rigid body alignment
2. **RMSD Calculation**: After alignment, RMSD is computed only on the designed protein chains

This approach ensures that we're measuring how consistently the designed binder maintains its binding pose relative to the target, rather than measuring the overall complex RMSD.

## Enabling Complex RMSD Metrics

To compute complex RMSD metrics during the analysis step, you need to enable them in the analysis configuration.

### Method 1: Using Command-Line Override

When running the full pipeline:

```bash
boltzgen run example/your_design.yaml \
  --output workbench/test_run \
  --protocol protein-anything \
  --num_designs 100 \
  --config analysis complex_rmsd_metrics=true complex_rmsd_threshold=2.5
```

### Method 2: Modifying Config File

Edit `config/analysis.yaml`:

```yaml
# Complex RMSD metrics (for filtering based on diffusion samples)
complex_rmsd_metrics: true  # Enable to compute complex RMSD across all diffusion samples
complex_rmsd_threshold: 2.5  # Threshold for counting passing diffusion samples
```

## Using Complex RMSD for Filtering

Once the complex RMSD metrics have been computed during analysis, you can filter designs based on them.

### Method 1: During Initial Pipeline Run

```bash
boltzgen run example/your_design.yaml \
  --output workbench/test_run \
  --protocol protein-anything \
  --num_designs 100 \
  --config analysis complex_rmsd_metrics=true complex_rmsd_threshold=2.5 \
  --config filtering filter_complex_rmsd=true min_complex_rmsd_samples=3
```

This will:
1. Compute complex RMSD for all diffusion samples during analysis
2. Filter to keep only designs where at least 3 diffusion samples pass the 2.5 Å threshold

### Method 2: Rerunning Only the Filtering Step

If you've already run the analysis with complex RMSD metrics enabled, you can quickly rerun just the filtering:

```bash
boltzgen run example/your_design.yaml \
  --output workbench/test_run \
  --protocol protein-anything \
  --steps filtering \
  --config filtering filter_complex_rmsd=true min_complex_rmsd_samples=5
```

### Method 3: Using the Jupyter Notebook

You can also use the provided `filter.ipynb` notebook for interactive filtering. After running analysis with complex RMSD metrics enabled:

1. Open `filter.ipynb`
2. Update the Filter initialization to include:
   ```python
   filter = Filter(
       design_dir="workbench/test_run",
       filter_complex_rmsd=True,
       min_complex_rmsd_samples=3,
       # ... other parameters
   )
   ```

## Configuration Parameters

### Analysis Parameters

- **`complex_rmsd_metrics`** (bool, default: `false`)
  - Enable computation of complex RMSD across all diffusion samples
  - Must be enabled to use complex RMSD filtering

- **`complex_rmsd_threshold`** (float, default: `2.5`)
  - RMSD threshold in Ångströms for counting passing samples
  - Samples with RMSD ≤ threshold are considered passing

### Filtering Parameters

- **`filter_complex_rmsd`** (bool, default: `false`)
  - Enable filtering based on complex RMSD metrics
  - Requires that `complex_rmsd_metrics=true` was used during analysis

- **`min_complex_rmsd_samples`** (int, default: `1`)
  - Minimum number of diffusion samples that must pass the threshold
  - Designs with fewer passing samples will be filtered out

## Output Metrics

When complex RMSD metrics are enabled, the following columns will be added to the analysis CSV:

- **`complex_rmsd_best`**: The best (minimum) RMSD across all diffusion samples
- **`complex_rmsd_samples_passing`**: Number of samples with RMSD ≤ threshold
- **`complex_rmsd_samples_total`**: Total number of diffusion samples
- **`complex_rmsd_pass_fraction`**: Fraction of samples passing (passing/total)

## Example Workflow

Here's a complete example workflow using complex RMSD filtering:

### Step 1: Generate Designs with Complex RMSD Analysis

```bash
boltzgen run example/vanilla_protein/1g13prot.yaml \
  --output workbench/protein_design \
  --protocol protein-anything \
  --num_designs 10000 \
  --budget 30 \
  --config analysis complex_rmsd_metrics=true complex_rmsd_threshold=2.0
```

This generates 10,000 designs and computes complex RMSD metrics during analysis.

### Step 2: Apply Strict Filtering

```bash
boltzgen run example/vanilla_protein/1g13prot.yaml \
  --output workbench/protein_design \
  --protocol protein-anything \
  --steps filtering \
  --config filtering filter_complex_rmsd=true min_complex_rmsd_samples=5 \
  --budget 30
```

This filters to keep only designs where at least 5 diffusion samples achieve RMSD ≤ 2.0 Å.

### Step 3: Examine Results

The filtered results will be in:
- `workbench/protein_design/final_ranked_designs/final_30_designs/` - Final selected structures
- `workbench/protein_design/final_ranked_designs/all_designs_metrics.csv` - All metrics including complex RMSD
- `workbench/protein_design/final_ranked_designs/results_overview.pdf` - Summary visualizations

## Interpreting Results

### What Does "Passing Samples" Mean?

A design with `complex_rmsd_samples_passing=5` and `complex_rmsd_samples_total=10` means:
- 10 diffusion samples were generated during refolding
- 5 of those samples achieved an RMSD ≤ threshold when:
  - Aligning the target structure
  - Computing RMSD on the designed binder only

This indicates the binding conformation was successfully recapitulated in 50% of the samples.

### Recommended Thresholds

Based on the use case:

- **Strict filtering** (high confidence in binding pose):
  - `complex_rmsd_threshold=2.0`
  - `min_complex_rmsd_samples=5` (or 50% of total samples)

- **Moderate filtering** (balanced):
  - `complex_rmsd_threshold=2.5`
  - `min_complex_rmsd_samples=3` (or 30% of total samples)

- **Lenient filtering** (exploratory):
  - `complex_rmsd_threshold=3.0`
  - `min_complex_rmsd_samples=2` (or 20% of total samples)

## Troubleshooting

### Error: Column 'complex_rmsd_samples_passing' not found

This means complex RMSD metrics were not computed during analysis. Solutions:
1. Rerun analysis with `complex_rmsd_metrics=true`
2. Or rerun the full pipeline with the analysis config enabled

### Warning: No designs pass the filter

If no designs pass the complex RMSD filter:
1. Check your threshold values - they may be too strict
2. Examine the distribution in `all_designs_metrics.csv`
3. Consider relaxing `min_complex_rmsd_samples` or `complex_rmsd_threshold`

### Running Only Analysis Step

If you want to add complex RMSD metrics to existing results:

```bash
boltzgen run example/your_design.yaml \
  --output workbench/test_run \
  --protocol protein-anything \
  --steps analysis \
  --config analysis complex_rmsd_metrics=true complex_rmsd_threshold=2.5
```

Note: This requires that the folding step was completed and the NPZ files exist.

## Performance Considerations

Computing complex RMSD metrics adds minimal overhead:
- Uses existing refolded structures (no additional folding needed)
- Computation time: ~1-2 seconds per design
- The metrics are only computed once during analysis

The filtering step itself is very fast (~30 seconds for 100k designs).

## Technical Details

For developers interested in the implementation:

- **Alignment function**: Uses `weighted_rigid_align()` from `boltzgen.model.loss.diffusion`
- **RMSD computation**: Uses `compute_subset_rmsd()` from `boltzgen.model.loss.validation`
- **Key function**: `compute_complex_rmsd_all_samples()` in `src/boltzgen/task/analyze/analyze_utils.py`

The implementation ensures:
- All diffusion samples are evaluated (not just the best)
- Proper masking of designed vs. target chains
- Efficient batch processing of samples
