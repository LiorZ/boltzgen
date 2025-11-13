# Quick Start: Complex RMSD Filtering

## TL;DR

Filter designs based on how many diffusion samples achieve consistent binding poses.

## One-Command Setup

```bash
boltzgen run your_design.yaml \
  --output workbench/output \
  --protocol protein-anything \
  --num_designs 10000 \
  --config analysis complex_rmsd_metrics=true complex_rmsd_threshold=2.5 \
  --config filtering filter_complex_rmsd=true min_complex_rmsd_samples=3
```

## What This Does

1. **Generates designs** as usual
2. **Computes complex RMSD** for each diffusion sample:
   - Aligns target (non-designed) chains
   - Measures RMSD on designed chains only
3. **Counts passing samples** (RMSD ≤ 2.5 Å)
4. **Filters designs** keeping only those with ≥ 3 passing samples

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `complex_rmsd_metrics` | `false` | Enable complex RMSD computation (analysis step) |
| `complex_rmsd_threshold` | `2.5` | RMSD threshold in Ångströms |
| `filter_complex_rmsd` | `false` | Enable filtering on complex RMSD (filtering step) |
| `min_complex_rmsd_samples` | `1` | Minimum passing samples required |

## Recommended Settings

### Strict (High Confidence)
```bash
--config analysis complex_rmsd_threshold=2.0 \
--config filtering min_complex_rmsd_samples=5
```

### Moderate (Balanced)
```bash
--config analysis complex_rmsd_threshold=2.5 \
--config filtering min_complex_rmsd_samples=3
```

### Lenient (Exploratory)
```bash
--config analysis complex_rmsd_threshold=3.0 \
--config filtering min_complex_rmsd_samples=2
```

## Rerun Filtering Only

Already have results? Just rerun filtering:

```bash
boltzgen run your_design.yaml \
  --output workbench/output \
  --steps filtering \
  --config filtering filter_complex_rmsd=true min_complex_rmsd_samples=5
```

⚠️ **Note**: Analysis must have been run with `complex_rmsd_metrics=true`

## Check Results

Look for these new columns in `all_designs_metrics.csv`:
- `complex_rmsd_best` - Best RMSD across all samples
- `complex_rmsd_samples_passing` - How many samples passed
- `complex_rmsd_pass_fraction` - Fraction passing (0.0 to 1.0)

## Full Documentation

See [COMPLEX_RMSD_FILTERING.md](COMPLEX_RMSD_FILTERING.md) for:
- Detailed explanations
- Multiple workflow examples
- Troubleshooting guide
- Technical implementation details
