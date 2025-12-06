# Membership Analysis - Quick Start Guide

## Overview

The `membership.ipynb` notebook computes cluster membership probabilities for **both M34 and M2** using:
- **CMD matching** (Color-Magnitude Diagram with isochrones)
- **Spatial distribution** (radial profiles)
- **Proper motion** (from Gaia, when available)

## Quick Start

### 1. Select Your Cluster

In **Cell 1**, set:
```python
CLUSTER = "M34"  # or "M2"
```

### 2. Restart Kernel & Run All Cells

The notebook will automatically:
- Load the correct isochrone
- Load all photometry from your .fits files
- Compute membership probabilities
- Save results to CSV

### 3. Output Files

The notebook saves:
- `m34_membership.csv` (for M34)
- `m2_membership.csv` (for M2)

Each file contains:
- All original photometry columns
- `P_cmd`: CMD membership probability
- `P_spatial`: Spatial membership probability
- `P_combined`: Combined membership probability (CMD × Spatial × PM if available)
- `has_pm`: Boolean flag for stars with Gaia proper motion data (if Gaia available)

## What Gets Loaded

| Cluster | Isochrone | Photometry | Stars | Age | Parameters |
|---------|-----------|------------|-------|-----|------------|
| **M34** | `isochrone_m34_corrected.dat` | `m34_photometry_calibrated.csv` | All from .fits | 0.2 Gyr | r_core=5.0' |
| **M2** | `isochrone_m2_13gyr_mh-1.6_corrected.dat` | `m2_photometry_calibrated.csv` | All from .fits | 13.0 Gyr | r_core=0.5' |

The photometry files contain **all stars extracted from your FITS images** during the data reduction pipeline.

## Working Without Gaia (If Gaia is Down)

The notebook **automatically detects** if Gaia data is unavailable and:
- ✓ Still computes CMD + Spatial membership
- ✗ Skips proper motion analysis
- ✓ Saves results with 2D membership (CMD × Spatial)

Output will show:
```
✗ Gaia data not available (Gaia may be down)
  Using CMD + Spatial only (no proper motion)
  Combination: CMD × Spatial only
```

## Interpreting Results

### Membership Probability Thresholds

| P_combined | Interpretation |
|------------|----------------|
| > 0.9 | High-confidence cluster member |
| 0.5 - 0.9 | Likely member |
| 0.1 - 0.5 | Uncertain |
| < 0.1 | Likely field star |

### Expected Results

**M34 (Open Cluster):**
- ~200-400 likely members (P > 0.5)
- ~50-150 high-confidence (P > 0.9)
- Most members at r < 10 arcmin

**M2 (Globular Cluster):**
- ~1000-3000 likely members (P > 0.5)
- ~500-1500 high-confidence (P > 0.9)
- Highly concentrated (r < 5 arcmin)

## Next Steps: Density Profiles

After getting membership results, you can:

1. **Extract members:**
   ```python
   members = pd.read_csv('m34_membership.csv')
   high_conf = members[members['P_combined'] > 0.9]
   ```

2. **Compute radial bins:**
   ```python
   r_bins = np.logspace(-1, 1.5, 20)  # 0.1 to 30 arcmin
   ```

3. **Surface density:**
   ```python
   N, _ = np.histogram(high_conf['r_arcmin'], bins=r_bins)
   area = np.pi * (r_bins[1:]**2 - r_bins[:-1]**2)
   sigma = N / area
   ```

4. **Fit Plummer profile:**
   ```python
   # σ(r) = σ_0 / (1 + (r/r_c)^2)^(5/2)
   ```

## Troubleshooting

### Issue: "NameError: name 'color_iso' is not defined"
**Solution**: You didn't run the isochrone loading cell. Run cells in order.

### Issue: "FileNotFoundError: m34_photometry_calibrated.csv"
**Solution**: Run `data_reduction_simple.ipynb` first to generate photometry catalogs.

### Issue: "Gaia query timeout"
**Solution**: Gaia archive may be down. The notebook will skip PM and use CMD+Spatial only.

### Issue: "No valid P_CMD values"
**Solution**: Check that your isochrone matches your data. Look at the CMD plot - does the isochrone overlay the data?

## Files Required

Before running, ensure you have:
- ✓ `m34_photometry_calibrated.csv` (from data reduction)
- ✓ `m2_photometry_calibrated.csv` (from data reduction)
- ✓ `isochrone_m34_corrected.dat` (PARSEC isochrone)
- ✓ `isochrone_m2_13gyr_mh-1.6_corrected.dat` (PARSEC isochrone)

All of these should already exist in your directory.

## Summary

This notebook:
1. Loads isochrone + photometry for selected cluster
2. Computes CMD membership (isochrone distance)
3. Computes spatial membership (Plummer profile)
4. (If Gaia up) Queries Gaia, cross-matches, adds PM membership
5. Combines all probabilities → P_combined
6. Saves results to CSV with all membership probabilities

You can now use the output CSV files for density profile analysis!
