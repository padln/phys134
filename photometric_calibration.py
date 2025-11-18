#!/usr/bin/env python3
"""
Photometric Calibration Script for M34

This script performs Steps 8.1-8.5 from data_reduction_simple.ipynb:
1. Query Gaia DR3 for M34 field
2. Convert Gaia photometry to SDSS
3. Cross-match with instrumental catalog
4. Compute zeropoints
5. Apply to all sources and save calibrated catalog

Run this if the notebook cells are giving you trouble with variable order.
"""

import numpy as np
import pandas as pd
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.gaia import Gaia
from astropy.stats import sigma_clip
import matplotlib.pyplot as plt

print("=" * 70)
print("M34 PHOTOMETRIC CALIBRATION")
print("=" * 70)

# ============================================================================
# Step 1: Load Instrumental Catalog
# ============================================================================
print("\n[1/5] Loading instrumental photometry catalog...")

try:
    catalog = Table.read('m34_photometry_instrumental.csv', format='csv')
    print(f"✓ Loaded {len(catalog)} sources")
    print(f"  Columns: {catalog.colnames}")
except FileNotFoundError:
    print("✗ ERROR: m34_photometry_instrumental.csv not found!")
    print("  Run data_reduction_simple.ipynb cells 1-16 first to create this file.")
    exit(1)

# ============================================================================
# Step 2: Query Gaia DR3
# ============================================================================
print("\n[2/5] Querying Gaia DR3 for M34 field...")

ra_M34 = 40.5092  # degrees
dec_M34 = 42.7619  # degrees
radius = 30.0 / 60.0  # 30 arcmin

query = f"""
SELECT source_id, ra, dec,
       phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag,
       pmra, pmdec, pmra_error, pmdec_error,
       parallax, parallax_error,
       phot_g_mean_flux_over_error,
       ruwe
FROM gaiadr3.gaia_source
WHERE 1=CONTAINS(
    POINT('ICRS', ra, dec),
    CIRCLE('ICRS', {ra_M34}, {dec_M34}, {radius})
)
AND phot_g_mean_mag IS NOT NULL
AND phot_bp_mean_mag IS NOT NULL
AND phot_rp_mean_mag IS NOT NULL
AND ruwe < 1.4
AND phot_g_mean_flux_over_error > 10
"""

job = Gaia.launch_job_async(query)
gaia_catalog = job.get_results()

print(f"✓ Retrieved {len(gaia_catalog)} Gaia sources")

# ============================================================================
# Step 3: Convert Gaia to SDSS
# ============================================================================
print("\n[3/5] Converting Gaia photometry to SDSS system...")

def gaia_to_sdss(G, BP, RP):
    """Convert Gaia G, BP, RP to SDSS g, r using Evans+ 2018 transformations."""
    BP_RP = BP - RP

    # g-band (Evans et al. 2018, Eq. 5.7)
    a0, a1, a2, a3 = -0.01760, -0.00686, -0.1069, -0.01256
    g_sdss = G + a0 + a1*BP_RP + a2*BP_RP**2 + a3*BP_RP**3

    # r-band (Evans et al. 2018, Eq. 5.8)
    b0, b1, b2, b3 = -0.03226, -0.3884, -0.1158, -0.01442
    r_sdss = G + b0 + b1*BP_RP + b2*BP_RP**2 + b3*BP_RP**3

    return g_sdss, r_sdss

gaia_g, gaia_r = gaia_to_sdss(
    gaia_catalog['phot_g_mean_mag'],
    gaia_catalog['phot_bp_mean_mag'],
    gaia_catalog['phot_rp_mean_mag']
)

gaia_catalog['g_sdss'] = gaia_g
gaia_catalog['r_sdss'] = gaia_r
gaia_catalog['g_r_sdss'] = gaia_g - gaia_r

print(f"✓ Converted {len(gaia_catalog)} sources to SDSS")
print(f"  SDSS g range: {np.min(gaia_g):.2f} - {np.max(gaia_g):.2f} mag")
print(f"  SDSS r range: {np.min(gaia_r):.2f} - {np.max(gaia_r):.2f} mag")

# ============================================================================
# Step 4: Cross-match and Compute Zeropoints
# ============================================================================
print("\n[4/5] Cross-matching catalogs and computing zeropoints...")

# Create SkyCoord objects
coords_inst = SkyCoord(ra=catalog['ra'], dec=catalog['dec'], unit='deg')
coords_gaia = SkyCoord(ra=gaia_catalog['ra'], dec=gaia_catalog['dec'], unit='deg')

# Match within 1 arcsec
idx_inst, idx_gaia, sep, _ = coords_gaia.search_around_sky(coords_inst, 1.0*u.arcsec)

print(f"✓ Matched {len(idx_inst)} sources (sep < 1 arcsec)")
print(f"  Median separation: {np.median(sep.arcsec):.3f} arcsec")

# Extract matched magnitudes
matched_gaia_g = gaia_catalog['g_sdss'][idx_gaia]
matched_gaia_r = gaia_catalog['r_sdss'][idx_gaia]
matched_inst_g = catalog['g_inst'][idx_inst]
matched_inst_r = catalog['r_inst'][idx_inst]
matched_g_err = catalog['g_err'][idx_inst]
matched_r_err = catalog['r_err'][idx_inst]
matched_color_gaia = gaia_catalog['g_r_sdss'][idx_gaia]

# Quality cuts for calibration stars
good_calib = (
    (matched_gaia_g > 12) & (matched_gaia_g < 18) &
    (matched_g_err < 0.05) & (matched_r_err < 0.05) &
    (matched_color_gaia > 0.2) & (matched_color_gaia < 1.5)
)

n_calib = np.sum(good_calib)
print(f"✓ Using {n_calib} stars for zeropoint (quality cuts applied)")

# Compute zeropoints with sigma clipping
zp_g_values = matched_gaia_g[good_calib] - matched_inst_g[good_calib]
zp_g_clipped = sigma_clip(zp_g_values, sigma=3, maxiters=5)
ZP_g = np.median(zp_g_clipped[~zp_g_clipped.mask])
ZP_g_std = np.std(zp_g_clipped[~zp_g_clipped.mask])

zp_r_values = matched_gaia_r[good_calib] - matched_inst_r[good_calib]
zp_r_clipped = sigma_clip(zp_r_values, sigma=3, maxiters=5)
ZP_r = np.median(zp_r_clipped[~zp_r_clipped.mask])
ZP_r_std = np.std(zp_r_clipped[~zp_r_clipped.mask])

print("\n" + "=" * 70)
print("PHOTOMETRIC ZEROPOINTS")
print("=" * 70)
print(f"g-band: ZP_g = {ZP_g:.3f} ± {ZP_g_std:.3f} mag")
print(f"r-band: ZP_r = {ZP_r:.3f} ± {ZP_r_std:.3f} mag")
print(f"Stars used: {np.sum(~zp_g_clipped.mask)} (g), {np.sum(~zp_r_clipped.mask)} (r)")
print("=" * 70)

# ============================================================================
# Step 5: Apply Calibration to All Sources
# ============================================================================
print("\n[5/5] Applying calibration to all sources...")

catalog['g'] = catalog['g_inst'] + ZP_g
catalog['r'] = catalog['r_inst'] + ZP_r
catalog['g_minus_r'] = catalog['g'] - catalog['r']

# Uncertainties (add zeropoint uncertainty in quadrature)
catalog['g_err_calib'] = np.sqrt(catalog['g_err']**2 + ZP_g_std**2)
catalog['r_err_calib'] = np.sqrt(catalog['r_err']**2 + ZP_r_std**2)

print(f"✓ Calibrated {len(catalog)} sources")
print(f"\nCalibrated magnitude range:")
print(f"  g: {catalog['g'].min():.2f} - {catalog['g'].max():.2f} mag")
print(f"  r: {catalog['r'].min():.2f} - {catalog['r'].max():.2f} mag")
print(f"  g-r: {catalog['g_minus_r'].min():.2f} - {catalog['g_minus_r'].max():.2f} mag")

# Save calibrated catalog
output_file = "m34_photometry_calibrated.csv"
catalog.write(output_file, format='csv', overwrite=True)
print(f"\n✓ Saved calibrated catalog to {output_file}")

# ============================================================================
# Diagnostic Plots
# ============================================================================
print("\n[OPTIONAL] Generating diagnostic plots...")

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Panel 1: Zeropoint residuals (g-band)
ax = axes[0, 0]
ax.scatter(matched_gaia_g[good_calib], zp_g_values, s=10, alpha=0.5, c='blue', label='All calib stars')
ax.scatter(matched_gaia_g[good_calib][~zp_g_clipped.mask],
           zp_g_values[~zp_g_clipped.mask], s=15, alpha=0.7, c='green', label='Used for ZP')
ax.axhline(ZP_g, color='red', linestyle='--', lw=2, label=f'ZP_g = {ZP_g:.3f}')
ax.axhline(ZP_g + ZP_g_std, color='red', linestyle=':', lw=1, alpha=0.5)
ax.axhline(ZP_g - ZP_g_std, color='red', linestyle=':', lw=1, alpha=0.5)
ax.set_xlabel('Gaia g [mag]')
ax.set_ylabel('g_Gaia - g_inst [mag]')
ax.set_title('g-band Zeropoint')
ax.legend()
ax.grid(alpha=0.3)

# Panel 2: Zeropoint residuals (r-band)
ax = axes[0, 1]
ax.scatter(matched_gaia_r[good_calib], zp_r_values, s=10, alpha=0.5, c='red', label='All calib stars')
ax.scatter(matched_gaia_r[good_calib][~zp_r_clipped.mask],
           zp_r_values[~zp_r_clipped.mask], s=15, alpha=0.7, c='green', label='Used for ZP')
ax.axhline(ZP_r, color='darkred', linestyle='--', lw=2, label=f'ZP_r = {ZP_r:.3f}')
ax.axhline(ZP_r + ZP_r_std, color='darkred', linestyle=':', lw=1, alpha=0.5)
ax.axhline(ZP_r - ZP_r_std, color='darkred', linestyle=':', lw=1, alpha=0.5)
ax.set_xlabel('Gaia r [mag]')
ax.set_ylabel('r_Gaia - r_inst [mag]')
ax.set_title('r-band Zeropoint')
ax.legend()
ax.grid(alpha=0.3)

# Panel 3: Instrumental CMD
ax = axes[1, 0]
ax.scatter(catalog['g_inst'] - catalog['r_inst'], catalog['g_inst'],
           s=5, alpha=0.5, c='gray')
ax.set_xlabel('g - r (instrumental)')
ax.set_ylabel('g (instrumental)')
ax.set_title('Instrumental CMD')
ax.invert_yaxis()
ax.grid(alpha=0.3)

# Panel 4: Calibrated CMD
ax = axes[1, 1]
ax.scatter(catalog['g_minus_r'], catalog['g'], s=5, alpha=0.5, c='black')
ax.set_xlabel('g - r [mag]')
ax.set_ylabel('g [mag]')
ax.set_title('Calibrated CMD (SDSS)')
ax.invert_yaxis()
ax.grid(alpha=0.3)

# Add zeropoint info
textstr = f'ZP_g = {ZP_g:.2f} ± {ZP_g_std:.2f}\nZP_r = {ZP_r:.2f} ± {ZP_r_std:.2f}'
ax.text(0.05, 0.95, textstr, transform=ax.transAxes,
        fontsize=10, verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.tight_layout()
plt.savefig('photometric_calibration_diagnostics.png', dpi=150, bbox_inches='tight')
print(f"✓ Saved diagnostic plots to photometric_calibration_diagnostics.png")

print("\n" + "=" * 70)
print("PHOTOMETRIC CALIBRATION COMPLETE!")
print("=" * 70)
print(f"\nNext steps:")
print(f"  1. Open membership.ipynb")
print(f"  2. Load m34_photometry_calibrated.csv")
print(f"  3. Run membership analysis")
print("=" * 70)
