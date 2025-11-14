# Data Reduction Pipeline - Implementation Guide

## Current Status

I've created **`data_reduction.ipynb`** with sections 1-2 complete (setup and image loading).

## What's Included So Far

✅ **Section 1: Setup** - All imports for photometry pipeline
✅ **Section 2: Load FITS** - Read images, extract metadata, visualize

## What Still Needs to be Added

The notebook needs the following sections to be complete. I can continue adding these:

### Section 3: Background Estimation and Subtraction
```python
# Use Background2D with sigma-clipped median
# Account for spatial variations
# Create background-subtracted images
```

### Section 4: Source Detection
```python
# DAOStarFinder or IRAFStarFinder
# Threshold = 3-5σ above background
# Filter cosmic rays and artifacts
```

### Section 5: Aperture Photometry
```python
# Optimal aperture radius (1.5-2× FWHM)
# Sky annulus for local background
# Aperture corrections
# Uncertainty propagation (Poisson + read noise + background)
```

### Section 6: PSF Modeling (Optional for M34, Critical for M2)
```python
# Extract bright, isolated stars
# Build empirical PSF (ePSF)
# PSF photometry for crowded regions
```

### Section 7: Astrometric Cross-Match with Gaia
```python
# Query Gaia DR3 in M34 field
# Match sources using WCS + nearest neighbor
# Add proper motions to catalog
```

### Section 8: Photometric Calibration
```python
# Use Gaia G, BP, RP → SDSS g,r transformation
# Or use Pan-STARRS photometry
# Determine zeropoints: m_inst + ZP = m_calibrated
# Account for airmass extinction
```

### Section 9: Multi-Band Catalog Creation
```python
# Cross-match g and r detections
# Create final catalog with:
#   - RA, Dec (from WCS)
#   - g, r magnitudes + uncertainties
#   - g-r color
#   - Gaia proper motions
#   - Quality flags
```

### Section 10: Quality Assessment
```python
# Magnitude histograms
# Color-magnitude diagram
# Spatial distribution
# Comparison with literature
```

### Section 11: Fisher Information Analysis
```python
# Compute Fisher information matrix for current observations
# Determine optimal additional exposure times
# Simulate S/N improvement for additional observations
# Recommendation: Take more data? Which band?
```

## Key Implementation Details

### Background Estimation Strategy
For M34 (moderately crowded):
- Use `Background2D` with box size ~100 pixels (74 arcsec)
- Sigma-clip to remove stars (3σ, 5 iterations)
- Median filter to handle spatial variations

### Source Detection Parameters
- **Threshold**: 3σ above local background (conservative)
- **FWHM estimate**: 2.3-2.8\" / 0.74\"/pix ≈ 3.1-3.8 pixels
- **Sharpness**: 0.2-1.0 (reject cosmic rays and galaxies)
- **Roundness**: -1 to 1 (reject cosmic rays)

### Aperture Photometry Setup
- **Aperture radius**: 1.5 × FWHM ≈ 5 pixels (captures ~90% of flux)
- **Sky annulus**: inner=7 pix, outer=10 pix
- **Aperture correction**: Measure from curve of growth on bright stars

### Photometric Uncertainties
Total uncertainty in magnitudes:

$$\\sigma_m = \\frac{2.5}{\\ln(10)} \\frac{\\sigma_\\text{flux}}{\\text{flux}}$$

where:

$$\\sigma_\\text{flux}^2 = \\frac{N_\\star}{g} + N_\\text{pix}\\left(\\frac{\\sigma_\\text{sky}^2}{g^2} + \\frac{R^2}{g^2}\\right)$$

- $N_\\star$: source counts (ADU)
- $N_\\text{pix}$: number of pixels in aperture
- $\\sigma_\\text{sky}$: sky background std (ADU)
- $R$: read noise (electrons)
- $g$: gain (e⁻/ADU)

### Gaia Cross-Matching
1. Query Gaia DR3 within 15' of M34 center (RA=02:42:09, Dec=+42:44:00)
2. Use astropy WCS to convert pixel → RA/Dec
3. Match with `match_coordinates_sky` (tolerance = 1-2 arcsec)
4. Keep only matches with RUWE < 1.4 (good astrometry)

### Photometric Calibration Options

**Option A: Gaia Photometry**
```python
# Transform Gaia G, BP, RP to SDSS g, r
# Use color transformations from Evans et al. (2018)
# Typical zeropoint accuracy: ~0.05 mag
```

**Option B: Pan-STARRS**
```python
# Query Pan-STARRS DR2 photometry
# Direct g, r magnitudes (native SDSS-like system)
# Better accuracy: ~0.02 mag
```

**Option C: Differential Photometry**
```python
# Use known M34 members from literature
# Fit zeropoint to match literature magnitudes
# Good for relative photometry
```

### Fisher Information for Observation Planning

The Fisher Information Matrix tells us how much \"information\" our observations contain about stellar parameters. For a single star with magnitude $m$ and error $\\sigma_m$:

$$I(m) = \\frac{1}{\\sigma_m^2}$$

For the entire cluster, we sum over all stars:

$$I_\\text{total} = \\sum_i \\frac{1}{\\sigma_{m,i}^2}$$

**Use case**:
- Compute $I_\\text{current}$ from existing data
- Simulate $I_\\text{new}$ with additional exposures
- Determine diminishing returns point

We can also compute Fisher information for parameters like:
- Plummer radius $a$
- Total mass $M$
- Distance modulus $\\mu$

## Running the Pipeline

### Expected Runtime
- Load images: <5 sec
- Background subtraction: ~10 sec
- Source detection: ~20 sec (per band)
- Aperture photometry: ~30 sec (per band)
- Gaia query + matching: ~30 sec
- Calibration: <5 sec
- Catalog generation: <5 sec

**Total: ~2-3 minutes**

### Expected Output

**Final catalog** (`m34_photometry.csv`):
```
   id      ra        dec       x     y    g_inst  g_err  r_inst  r_err    g      r     g-r   pmra  pmdec  pmra_err pmdec_err  ...
0  1    40.5123   42.7334   1200  1150  18.234  0.015  17.891  0.012  18.24  17.89  0.35  -5.32  -3.41   0.12      0.15  ...
1  2    40.5089   42.7389   1205  1180  19.567  0.031  19.102  0.025  19.57  19.10  0.47  -5.18  -3.55   0.21      0.24  ...
...
```

~5000-8000 sources detected, of which ~200-400 are likely M34 members.

## Next Steps After Pipeline Completion

1. **Run artificial star tests** using `completeness.ipynb`
2. **Apply membership determination** using `membership_enhanced.ipynb`
3. **Construct density profiles** using new `density_profile.ipynb`
4. **Fit Plummer model** and compare to M2

## How to Continue Development

I can add the remaining sections (3-11) to the notebook. Would you like me to:

**Option A**: Add all sections at once (complete pipeline, ~1500 lines of code)

**Option B**: Add sections incrementally so you can test each step (recommended)

**Option C**: Add specific sections you need most urgently (e.g., just photometry + catalog)

Let me know and I'll continue building!
