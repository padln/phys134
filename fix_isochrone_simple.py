#!/usr/bin/env python3
"""
Simple Isochrone Correction for M34

Problem: PARSEC isochrone is at 10 pc with E(B-V) ~ 0.93
Solution: Move to 470 pc and apply E(B-V) = 0.08

This version is simpler and more robust.
"""

import numpy as np
import pandas as pd

print("=" * 70)
print("ISOCHRONE CORRECTION FOR M34")
print("=" * 70)

# Parameters
E_BV_remove = 0.93  # Current (wrong) reddening in file
E_BV_M34 = 0.08     # Correct M34 reddening
distance_M34 = 470  # M34 distance in pc
mu_M34 = 5 * np.log10(distance_M34 / 10)  # Distance modulus from 10 pc

# Extinction coefficients (Cardelli+ 1989, R_V=3.1)
A_g_per_EBV = 3.303
A_r_per_EBV = 2.285

print(f"\nInput file: isochrone_m34.dat (10 pc, E(B-V) ≈ {E_BV_remove})")
print(f"Output file: isochrone_m34_corrected.dat ({distance_M34} pc, E(B-V) = {E_BV_M34})")

# Load original isochrone
column_names = [
    'Zini', 'MH', 'logAge', 'Mini', 'int_IMF', 'Mass', 'logL', 'logTe',
    'logg', 'label', 'McoreTP', 'C_O', 'period0', 'period1', 'period2',
    'period3', 'period4', 'pmode', 'Mloss', 'tau1m', 'X', 'Y', 'Xc',
    'Xn', 'Xo', 'Cexcess', 'Z', 'mbolmag', 'umag', 'gmag', 'rmag',
    'imag', 'zmag'
]

iso = pd.read_csv('isochrone_m34.dat', delim_whitespace=True, comment='#',
                  names=column_names, na_values=['-1.00', '-9.999'])

print(f"✓ Loaded {len(iso)} isochrone points")

# Extract magnitudes as numpy arrays (avoid pandas weirdness)
g_orig = iso['gmag'].values.astype(float)
r_orig = iso['rmag'].values.astype(float)
i_orig = iso['imag'].values.astype(float)
z_orig = iso['zmag'].values.astype(float)
u_orig = iso['umag'].values.astype(float)

print(f"\nOriginal magnitudes (10 pc, E(B-V)={E_BV_remove}):")
print(f"  g: {g_orig.min():.2f} - {g_orig.max():.2f}")
print(f"  (g-r): {(g_orig - r_orig).min():.2f} - {(g_orig - r_orig).max():.2f}")

# Step 1: Remove current reddening
A_g_remove = A_g_per_EBV * E_BV_remove
A_r_remove = A_r_per_EBV * E_BV_remove
A_i_remove = 1.698 * E_BV_remove  # SDSS i-band
A_z_remove = 1.263 * E_BV_remove  # SDSS z-band
A_u_remove = 4.239 * E_BV_remove  # SDSS u-band

g_intrinsic_10pc = g_orig - A_g_remove
r_intrinsic_10pc = r_orig - A_r_remove
i_intrinsic_10pc = i_orig - A_i_remove
z_intrinsic_10pc = z_orig - A_z_remove
u_intrinsic_10pc = u_orig - A_u_remove

print(f"\nAfter removing reddening (10 pc, E(B-V)=0):")
print(f"  g: {g_intrinsic_10pc.min():.2f} - {g_intrinsic_10pc.max():.2f}")
print(f"  (g-r): {(g_intrinsic_10pc - r_intrinsic_10pc).min():.2f} - {(g_intrinsic_10pc - r_intrinsic_10pc).max():.2f}")

# Step 2: Apply distance modulus
g_intrinsic_M34 = g_intrinsic_10pc + mu_M34
r_intrinsic_M34 = r_intrinsic_10pc + mu_M34
i_intrinsic_M34 = i_intrinsic_10pc + mu_M34
z_intrinsic_M34 = z_intrinsic_10pc + mu_M34
u_intrinsic_M34 = u_intrinsic_10pc + mu_M34

print(f"\nAfter distance correction ({distance_M34} pc, E(B-V)=0):")
print(f"  Distance modulus: μ = {mu_M34:.3f} mag")
print(f"  g: {g_intrinsic_M34.min():.2f} - {g_intrinsic_M34.max():.2f}")

# Step 3: Apply M34 reddening
A_g_M34 = A_g_per_EBV * E_BV_M34
A_r_M34 = A_r_per_EBV * E_BV_M34
A_i_M34 = 1.698 * E_BV_M34
A_z_M34 = 1.263 * E_BV_M34
A_u_M34 = 4.239 * E_BV_M34

g_final = g_intrinsic_M34 + A_g_M34
r_final = r_intrinsic_M34 + A_r_M34
i_final = i_intrinsic_M34 + A_i_M34
z_final = z_intrinsic_M34 + A_z_M34
u_final = u_intrinsic_M34 + A_u_M34

print(f"\nFinal ({distance_M34} pc, E(B-V)={E_BV_M34}):")
print(f"  g: {g_final.min():.2f} - {g_final.max():.2f}")
print(f"  (g-r): {(g_final - r_final).min():.2f} - {(g_final - r_final).max():.2f}")

# Create corrected isochrone (make a clean copy, don't modify in place)
iso_corrected = iso.copy()

# Replace magnitudes with corrected values (force to numpy arrays to avoid pandas issues)
iso_corrected['gmag'] = g_final
iso_corrected['rmag'] = r_final
iso_corrected['imag'] = i_final
iso_corrected['zmag'] = z_final
iso_corrected['umag'] = u_final

# Save with space separation (like original)
output_file = 'isochrone_m34_corrected.dat'
iso_corrected.to_csv(output_file, sep=' ', index=False, float_format='%.6f')

print(f"\n✓ Saved corrected isochrone to {output_file}")

# Verification
print("\n" + "=" * 70)
print("VERIFICATION")
print("=" * 70)

# Check a solar-type star (0.5 M_sun)
mask_solar = (iso_corrected['Mini'] > 0.45) & (iso_corrected['Mini'] < 0.55) & (iso_corrected['label'] == 0)
if np.any(mask_solar):
    idx = np.where(mask_solar)[0][0]
    g_sol = iso_corrected['gmag'].iloc[idx]
    r_sol = iso_corrected['rmag'].iloc[idx]
    color_sol = g_sol - r_sol
    print(f"Solar-type star (0.5 M_sun):")
    print(f"  g = {g_sol:.2f} mag")
    print(f"  (g-r) = {color_sol:.2f} mag")
    print(f"  Expected: (g-r) ~ 0.5-0.6 mag for slightly reddened solar-type")
    if 0.4 < color_sol < 0.7:
        print(f"  ✓ Color is correct!")
    else:
        print(f"  ✗ Color seems wrong")

print("\n" + "=" * 70)
print("COMPLETE!")
print("=" * 70)
