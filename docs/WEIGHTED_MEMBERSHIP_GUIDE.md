# Weighted Membership Approach - No Hard Cuts!

## Philosophy

Instead of using arbitrary probability thresholds (P > 0.5 or P > 0.9), we use **continuous membership probabilities** as statistical weights. This:

1. **Preserves information** - Every star contributes according to its likelihood of membership
2. **Avoids arbitrary cutoffs** - No need to decide "is 0.49 a member or not?"
3. **Reduces bias** - Stars near the threshold don't cause discontinuities
4. **Proper statistics** - Weighted counts naturally incorporate uncertainty

## Key Changes Made

### 1. Fixed Binary Probability Issue

**Before** (WRONG):
```python
P_combined = P_combined / (P_combined + 1e-10)
```
This collapsed everything to 0 or 1!

**After** (CORRECT):
```python
P_combined = P_cmd * P_spatial  # Simple multiplication, no normalization!
```
Keeps continuous probabilities in [0,1] range.

### 2. Improved Spatial Background Estimation

**Before**: Fixed 10% of peak
**After**: Estimated from outer regions (r > 3 × r_core)

This gives better contrast and shows central concentration more clearly.

### 3. Weighted Counting

**Instead of**:
```python
n_members = np.sum(P_combined > 0.5)  # Hard cut
```

**Use**:
```python
weighted_members = np.sum(P_combined)  # Continuous weights
```

## How to Use Weighted Membership

### For Total Member Count

```python
# Each star contributes its probability
total_members = np.sum(P_combined)

# Example: If you have 3 stars with P = [0.9, 0.7, 0.2]
# total_members = 0.9 + 0.7 + 0.2 = 1.8 stars
```

### For Radial Density Profiles

```python
# Bin stars by radius
r_bins = np.logspace(-1, 1.5, 20)
r_centers = np.sqrt(r_bins[1:] * r_bins[:-1])

# Count weighted by membership probability
N_weighted, _ = np.histogram(r_arcmin, bins=r_bins, weights=P_combined)

# Surface density
area = np.pi * (r_bins[1:]**2 - r_bins[:-1]**2)
sigma = N_weighted / area

# Now sigma(r) represents the true member density!
```

### For Mass Estimation

```python
# If you have stellar masses
masses = estimate_masses(stars_g)  # From mass-luminosity relation

# Total cluster mass
M_cluster = np.sum(masses * P_combined)

# Each star contributes: (its mass) × (probability it's a member)
```

### For Luminosity Functions

```python
# Magnitude bins
mag_bins = np.arange(12, 22, 0.5)

# Weighted histogram
LF, _ = np.histogram(stars_g, bins=mag_bins, weights=P_combined)

# This is the member luminosity function, naturally correcting for field contamination
```

## Interpreting P_combined Values

| P_combined | Meaning | How to Think About It |
|------------|---------|----------------------|
| 0.95 | 95% certain it's a member | Include with high confidence |
| 0.70 | Likely member, some uncertainty | Count as 0.7 of a star |
| 0.50 | Equal probability member/field | Ambiguous - contributes 50% |
| 0.20 | Likely field star | Count as 0.2 of a star |
| 0.05 | Almost certainly field | Basically ignore |

## Example: Comparing Hard Cut vs Weighted

Suppose you have 5 stars with P_combined = [0.95, 0.85, 0.75, 0.55, 0.45]

**Hard cut at P > 0.5:**
- Members: 4 stars (includes 0.55, excludes 0.45)
- Problem: Very similar stars (0.55 vs 0.45) treated completely differently!

**Weighted approach:**
- Total membership: 0.95 + 0.85 + 0.75 + 0.55 + 0.45 = 3.55 stars
- Benefit: All stars contribute fairly, uncertainty preserved

## Density Profile Example

```python
import numpy as np
import matplotlib.pyplot as plt

# Load membership results
import pandas as pd
members = pd.read_csv('m34_membership.csv')

# Extract data
r = members['r_arcmin'].values  # Radial distance from center
P = members['P_combined'].values  # Membership probability

# Radial bins (logarithmic for wide dynamic range)
r_bins = np.logspace(-1, 1.5, 20)  # 0.1 to 30 arcmin
r_centers = np.sqrt(r_bins[1:] * r_bins[:-1])

# Weighted histogram
N_weighted, _ = np.histogram(r, bins=r_bins, weights=P)

# Surface density
area = np.pi * (r_bins[1:]**2 - r_bins[:-1]**2)
sigma = N_weighted / area

# Plot
plt.figure(figsize=(10, 6))
plt.loglog(r_centers, sigma, 'o-', lw=2, label='M34 density profile')
plt.xlabel('Radius [arcmin]')
plt.ylabel(r'Surface density [stars/arcmin$^2$]')
plt.title('M34 Radial Density Profile (Weighted by Membership)')
plt.grid(True, alpha=0.3)
plt.legend()
plt.show()

# Fit Plummer profile
from scipy.optimize import curve_fit

def plummer(r, sigma0, rc):
    return sigma0 / (1 + (r/rc)**2)**(5/2)

# Fit (use only points with enough stars)
good = N_weighted > 5
popt, pcov = curve_fit(plummer, r_centers[good], sigma[good])

print(f"Best-fit Plummer profile:")
print(f"  σ_0 = {popt[0]:.1f} stars/arcmin²")
print(f"  r_c = {popt[1]:.2f} arcmin")
```

## Why This Works

The weighted approach is mathematically equivalent to:
1. Running many Monte Carlo realizations where each star is included/excluded based on its probability
2. Averaging the results over all realizations
3. But computationally much faster!

## Statistical Uncertainties

The uncertainty in weighted counts follows Poisson statistics:

```python
# Uncertainty in total member count
sigma_N = np.sqrt(np.sum(P_combined**2))

# This accounts for the fact that high-P stars contribute more
# to both the count AND its uncertainty
```

## When to Use Hard Cuts (Rarely!)

Hard cuts are only needed when:
- Creating a clean sample for spectroscopic follow-up (use P > 0.9)
- Making a publication-quality CMD of "members" (use P > 0.7 or 0.8)
- Visual inspection (highlight P > 0.5 stars)

For all quantitative analysis (densities, masses, luminosities), **use weighted probabilities**.

## Summary

✓ Keep continuous probabilities [0,1]
✓ Use P_combined as statistical weights
✓ No arbitrary thresholds needed
✗ Don't collapse to binary (member/non-member)
✗ Don't normalize probabilities incorrectly

This approach is standard in modern cluster analysis (e.g., Gaia Collaboration papers on open clusters).
