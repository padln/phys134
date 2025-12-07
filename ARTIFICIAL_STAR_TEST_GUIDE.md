# Artificial Star Test Implementation Guide

## Overview

This guide describes how to use the artificial star test (AST) framework to measure photometric completeness in M2 and M34 observations.

## What is an Artificial Star Test?

An artificial star test measures how efficiently your photometry pipeline detects stars at different brightness levels. This is critical for:
- Correcting observed star counts to true populations
- Understanding survey depth and detection limits
- Quantifying systematic biases in photometric measurements

## Quick Start

### Run the Example Notebook

```bash
jupyter notebook artificial_star_test_example.ipynb
```

This notebook will:
1. Load FITS images from the `Data/` directory
2. Inject artificial stars at known positions and magnitudes (instrumental range: -16 to -8 mag)
3. Run the photometry pipeline to detect sources
4. Calculate recovery fractions (completeness)
5. Fit a completeness model to the data
6. Visualize and save results

**Note**: The magnitude range is set to match the instrumental magnitudes of real stars in the images.
Typical bright stars have instrumental magnitudes around -15 mag, while fainter detectable stars
are around -8 mag (before photometric calibration).

### Using the Python API

```python
from pipeline_utils import run_artificial_star_test, load_fits_image
import numpy as np

# Load your FITS image
data, header, wcs = load_fits_image("Data/your_image.fits")

# Define magnitude bins to test
mag_bins = np.arange(18, 24, 0.5)  # Test from mag 18 to 24

# Run the artificial star test
results = run_artificial_star_test(
    data=data,
    header=header,
    wcs=wcs,
    magnitude_bins=mag_bins,
    n_stars_per_bin=100,  # 100 stars per magnitude
    fwhm=3.5,  # PSF FWHM in pixels
    zeropoint=0.0,  # Use instrumental mags
    match_radius=2.0  # Match within 2 pixels
)

# View results
from pipeline_utils import print_ast_results
print_ast_results(results)
```

## Core Functions

### `add_artificial_stars(data, positions, magnitudes, fwhm, zeropoint=0.0)`

Injects artificial stars into an image using a Gaussian PSF model.

**Parameters:**
- `data`: 2D image array
- `positions`: (N, 2) array of (x, y) pixel coordinates
- `magnitudes`: Array of magnitudes to inject
- `fwhm`: PSF FWHM in pixels
- `zeropoint`: Photometric zeropoint (default=0 for instrumental mags)

**Returns:**
- Copy of image with artificial stars added

### `run_artificial_star_test(data, header, wcs, magnitude_bins, n_stars_per_bin, fwhm, ...)`

Performs complete artificial star test workflow.

**Parameters:**
- `data`: Original image data
- `header`: FITS header
- `wcs`: World coordinate system
- `magnitude_bins`: Array of magnitudes to test
- `n_stars_per_bin`: Number of stars to inject per magnitude
- `fwhm`: PSF FWHM in pixels
- `zeropoint`: Photometric zeropoint (optional)
- `detection_params`: Dict of detection parameters (optional)
- `photometry_params`: Dict of photometry parameters (optional)
- `match_radius`: Matching radius in pixels (default=3.0)

**Returns:**
Dictionary with:
- `magnitude_bins`: Input magnitude bins
- `n_added`: Stars added per bin
- `n_recovered`: Stars recovered per bin
- `completeness`: Recovery fraction per bin
- `completeness_err`: Binomial uncertainty per bin
- `injected_positions`: All injected positions
- `injected_magnitudes`: All injected magnitudes

### `print_ast_results(results)`

Prints formatted table of AST results.

## Completeness Model Fitting

After running the AST, fit a completeness model to the data:

```python
from scipy.special import erf
from scipy.optimize import minimize
import numpy as np

# Define error function model
def completeness_model(m, m50, sigma_comp):
    """Error function completeness model"""
    return 0.5 * (1 + erf((m50 - m) / (np.sqrt(2) * sigma_comp)))

# Negative log-likelihood
def negative_log_likelihood(theta, m_bins, N_add, N_rec):
    m50, sigma_comp = theta
    C = completeness_model(m_bins, m50, sigma_comp)
    C = np.clip(C, 1e-10, 1 - 1e-10)
    log_L = np.sum(N_rec * np.log(C) + (N_add - N_rec) * np.log(1 - C))
    return -log_L

# Fit the model
result_fit = minimize(
    negative_log_likelihood,
    [21.0, 1.0],  # Initial guess
    args=(results['magnitude_bins'], 
          results['n_added'], 
          results['n_recovered']),
    method='Nelder-Mead'
)

m50, sigma = result_fit.x
print(f"50% completeness at magnitude: {m50:.2f}")
print(f"Completeness width: {sigma:.2f}")
```

## Interpreting Results

### Key Parameters

- **m50**: Magnitude at 50% completeness - the limiting magnitude of your survey
- **sigma_comp**: Controls steepness of completeness transition
  - Smaller σ → sharper cutoff (better detector/seeing)
  - Larger σ → gradual decline (crowding, variable seeing)

### Completeness Levels

- **90% completeness**: m₅₀ + 1.28σ
- **50% completeness**: m₅₀
- **10% completeness**: m₅₀ - 1.28σ

### Using Completeness for Corrections

Correct observed star counts to true populations:

```python
N_true(m) = N_observed(m) / C(m)
```

Weight cluster membership probabilities:

```python
P_member_weighted = P_member(CMD, PM, spatial) / C(m)
```

### Using the Completeness Correction Module

```python
from completeness_correction import (
    load_completeness_model,
    apply_completeness_correction
)

# Load AST results
m50, sigma_comp = load_completeness_model('completeness_model_params.txt')

# Apply correction to radial profile
N_corrected, C = apply_completeness_correction(
    N_observed, magnitudes, m50, sigma_comp
)
```

See `completeness.ipynb` for detailed examples and `membership.ipynb` for integration
with radial profile analysis.

## Best Practices

### Number of Artificial Stars

- **Minimum**: 20-30 stars per magnitude bin
- **Recommended**: 50-100 stars per bin
- **High precision**: 200+ stars per bin

More stars = better statistics but longer computation time.

### Magnitude Range

- Start ~2 magnitudes brighter than brightest real sources
- End at the faint detection limit
- Use 0.5-1.0 mag bins for smooth coverage
- **For instrumental magnitudes**: Typical range is -16 to -8 mag
- **For calibrated magnitudes**: Adjust based on photometric zeropoint

### Spatial Considerations

- Avoid injecting stars near image edges (50-100 pixel border)
- For crowded fields, consider spatial variations in completeness
- Test different regions separately if necessary

### PSF Modeling

- Use measured FWHM from actual sources in the image
- For better accuracy, consider using actual PSF models instead of Gaussian
- Account for PSF variations across the field

## Example Results

Testing with M2 r-band data (79s exposure):
- **350 artificial stars injected** across 7 magnitude bins
- **179 stars recovered** (51% overall)
- **50% completeness at** instrumental magnitude -5.41
- **Completeness width** σ = 1.21 mag

Completeness curve shows expected behavior:
- Bright stars (mag < -6): >90% recovery
- Transition regime (mag -6 to -5): 50-90% recovery
- Faint stars (mag > -4): <20% recovery

## Troubleshooting

### Low recovery rates everywhere
- Check magnitude scale (instrumental vs calibrated)
- Verify FWHM is appropriate for image
- Adjust detection threshold parameters

### Unrealistic completeness (>100% or <0%)
- Usually a magnitude scale issue
- Check zeropoint parameter
- Verify flux-to-magnitude conversion

### Noisy completeness curve
- Increase n_stars_per_bin
- Use finer magnitude bins
- Check for crowding or artifacts

## References

- Stetson, P. B. 1994, PASP, 106, 250 (DAOPHOT photometry)
- Aparicio & Gallart 2004, AJ, 128, 1465 (Completeness corrections)
- Dolphin 2002, MNRAS, 332, 91 (Artificial star tests for HST)

## Citation

If you use this code in your research, please cite:

```bibtex
@misc{phys134_ast,
  author = {Nathan Madsen},
  title = {Artificial Star Test Framework for Stellar Cluster Photometry},
  year = {2025},
  publisher = {GitHub},
  url = {https://github.com/padln/phys134}
}
```

## Contact

For questions or issues:
- Open an issue on GitHub: https://github.com/padln/phys134
- Email: madsen@ucsb.edu

---

**Course**: PHYS 134L - Advanced Observational Astrophysics  
**Institution**: UC Santa Barbara  
**Term**: Fall 2025
