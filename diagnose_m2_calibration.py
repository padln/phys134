#!/usr/bin/env python3
"""
Diagnose M2 photometry calibration issue
"""

import pandas as pd
import numpy as np

# Load the calibrated catalog
catalog = pd.read_csv('m2_photometry_calibrated.csv')

print("=" * 70)
print("M2 PHOTOMETRY DIAGNOSIS")
print("=" * 70)

# Check magnitude ranges
print(f"\nCalibrated magnitude ranges:")
print(f"  g: {catalog['g'].min():.2f} to {catalog['g'].max():.2f}")
print(f"  r: {catalog['r'].min():.2f} to {catalog['r'].max():.2f}")
print(f"  g-r: {catalog['g_minus_r'].min():.2f} to {catalog['g_minus_r'].max():.2f}")

print(f"\nInstrumental magnitude ranges:")
print(f"  g_inst: {catalog['g_inst'].min():.2f} to {catalog['g_inst'].max():.2f}")
print(f"  r_inst: {catalog['r_inst'].min():.2f} to {catalog['r_inst'].max():.2f}")

# Infer zeropoint
brightest_idx = np.argmin(catalog['g'])
zp_g_inferred = catalog['g'].iloc[brightest_idx] - catalog['g_inst'].iloc[brightest_idx]
zp_r_inferred = catalog['r'].iloc[brightest_idx] - catalog['r_inst'].iloc[brightest_idx]

print(f"\nInferred zeropoints (from brightest star):")
print(f"  ZP_g = {zp_g_inferred:.2f} mag")
print(f"  ZP_r = {zp_r_inferred:.2f} mag")

print(f"\n--- PROBLEM IDENTIFIED ---")
print(f"For M2 at 11,500 pc (μ = 15.30 mag):")
print(f"  Expected RGB tip: g ~ 13-14 mag")
print(f"  Expected HB stars: g ~ 15-16 mag")
print(f"  Expected turnoff: g ~ 18-19 mag")
print(f"\nActual brightest stars: g ~ {catalog['g'].min():.1f} mag")
print(f"  → 3-4 magnitudes TOO BRIGHT!")

print(f"\n--- LIKELY CAUSE ---")
print(f"The Gaia calibration stars are FOREGROUND FIELD STARS, not M2 members.")
print(f"M2 is a distant globular cluster (11.5 kpc), so most Gaia matches")
print(f"will be nearby field stars at g~10-15, not cluster RGB at g~13-19.")

print(f"\n--- SOLUTION ---")
print(f"Recalibrate using:")
print(f"  1. Known M2 parameters (distance, reddening)")
print(f"  2. Standard star fields (if available)")
print(f"  3. Pan-STARRS catalog (better faint star coverage)")
print(f"  4. OR accept instrumental mags and use isochrone at same scale")

print("=" * 70)
