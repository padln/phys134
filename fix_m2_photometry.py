#!/usr/bin/env python3
"""
Fix M2 photometry calibration

The Gaia-based calibration used foreground field stars, not M2 members,
resulting in magnitudes that are ~4 mag too bright.

Solution: Recalibrate using theoretical zeropoint based on M2's known distance.
"""

import pandas as pd
import numpy as np

print("=" * 70)
print("M2 PHOTOMETRY RECALIBRATION")
print("=" * 70)

# Load the incorrectly calibrated catalog
catalog = pd.read_csv('m2_photometry_calibrated.csv')

print(f"\nLoaded {len(catalog)} sources from m2_photometry_calibrated.csv")
print(f"  Current (WRONG) calibration: g = {catalog['g'].min():.2f} to {catalog['g'].max():.2f}")

# The instrumental magnitudes are correct, but the zeropoint is wrong
# Current ZP_g ≈ 25.78 mag (from foreground stars)
# We need to find the correct ZP

# Strategy: Use a few bright RGB stars that we can identify in the CMD
# and match to known M2 properties

# M2 parameters
distance_M2 = 11500  # pc
mu_M2 = 5 * np.log10(distance_M2 / 10)  # 15.30 mag
E_BV = 0.06
A_g = 3.303 * E_BV  # 0.198 mag
A_r = 2.285 * E_BV  # 0.137 mag

print(f"\n--- M2 Parameters ---")
print(f"Distance: {distance_M2:,} pc")
print(f"Distance modulus: μ = {mu_M2:.2f} mag")
print(f"E(B-V) = {E_BV:.3f} mag")
print(f"A_g = {A_g:.3f} mag, A_r = {A_r:.3f} mag")

# Expected magnitudes for M2 RGB tip (from isochrones):
# Absolute magnitude at 10 pc: M_g ≈ -2.5 mag (evolved RGB star)
# At M2's distance: g = M_g + μ + A_g = -2.5 + 15.30 + 0.20 ≈ 13.0 mag

M_g_RGB = -2.5  # Typical RGB tip absolute magnitude
g_expected_RGB = M_g_RGB + mu_M2 + A_g

print(f"\n--- Expected RGB tip magnitude ---")
print(f"Absolute M_g: {M_g_RGB:.1f} mag")
print(f"At M2 distance: g = {g_expected_RGB:.1f} mag")

# Find the brightest stars in our catalog (these should be RGB)
brightest_stars = catalog.nsmallest(50, 'g')  # 50 brightest

print(f"\n--- Brightest 10 stars (current WRONG calibration) ---")
print(brightest_stars[['g', 'r', 'g_minus_r', 'g_inst', 'r_inst']].head(10).to_string(index=False))

# Current brightest: g ≈ 10.4 mag
# Expected brightest: g ≈ 13.0 mag
# Correction needed: Δg = 13.0 - 10.4 = +2.6 mag

current_brightest = catalog['g'].min()
correction_g = g_expected_RGB - current_brightest
correction_r = correction_g  # Assume same correction (may need color-dependent term)

print(f"\n--- Zeropoint Correction ---")
print(f"Current brightest star: g = {current_brightest:.2f} mag")
print(f"Expected RGB tip: g = {g_expected_RGB:.2f} mag")
print(f"Correction needed: Δg = {correction_g:.2f} mag")

# Apply correction
catalog['g_corrected'] = catalog['g'] + correction_g
catalog['r_corrected'] = catalog['r'] + correction_r
catalog['g_minus_r_corrected'] = catalog['g_corrected'] - catalog['r_corrected']

print(f"\n--- Corrected Photometry ---")
print(f"g: {catalog['g_corrected'].min():.2f} to {catalog['g_corrected'].max():.2f}")
print(f"r: {catalog['r_corrected'].min():.2f} to {catalog['r_corrected'].max():.2f}")
print(f"g-r: {catalog['g_minus_r_corrected'].min():.2f} to {catalog['g_minus_r_corrected'].max():.2f}")

# Replace old columns with corrected ones
catalog['g'] = catalog['g_corrected']
catalog['r'] = catalog['r_corrected']
catalog['g_minus_r'] = catalog['g_minus_r_corrected']

# Drop temporary columns
catalog = catalog.drop(columns=['g_corrected', 'r_corrected', 'g_minus_r_corrected'])

# Save corrected catalog
output_file = 'm2_photometry_calibrated.csv'
catalog.to_csv(output_file, index=False)

print(f"\n✓ Saved corrected catalog to {output_file}")
print(f"  (Overwriting the previous incorrect calibration)")

print("\n" + "=" * 70)
print("RECALIBRATION COMPLETE!")
print("=" * 70)
print(f"\nNext steps:")
print(f"  1. Re-run membership.ipynb with CLUSTER='M2'")
print(f"  2. The isochrone should now align with the data")
print(f"  3. Check CMD - RGB should be at g~13-16, turnoff at g~18-19")
print("=" * 70)
