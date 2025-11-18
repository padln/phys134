#!/usr/bin/env python3
"""
De-redden PARSEC Isochrone and Re-apply Correct M34 Extinction

Problem: The downloaded isochrone has too much reddening applied.
Solution: Remove existing reddening, apply correct E(B-V) = 0.08
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

print("=" * 70)
print("ISOCHRONE REDDENING CORRECTION")
print("=" * 70)

# Load isochrone
print("\nLoading PARSEC isochrone...")
column_names = [
    'Zini', 'MH', 'logAge', 'Mini', 'int_IMF', 'Mass', 'logL', 'logTe',
    'logg', 'label', 'McoreTP', 'C_O', 'period0', 'period1', 'period2',
    'period3', 'period4', 'pmode', 'Mloss', 'tau1m', 'X', 'Y', 'Xc',
    'Xn', 'Xo', 'Cexcess', 'Z', 'mbolmag', 'umag', 'gmag', 'rmag',
    'imag', 'zmag'
]

iso = pd.read_csv('isochrone_m34.dat', delim_whitespace=True, comment='#',
                  names=column_names, na_values=['-1.00', '-9.999'])

# Filter to 200 Myr, [M/H] = -0.1
unique_ages = iso['logAge'].unique()
log_age_200 = unique_ages[np.argmin(np.abs(10**(unique_ages-6) - 200))]
iso_200 = iso[(iso['logAge'] == log_age_200) & (iso['MH'] == -0.1)].copy()

print(f"✓ Loaded {len(iso_200)} points for 200 Myr, [M/H]=-0.1")

# Current magnitudes (with unknown reddening)
g_current = iso_200['gmag'].values
r_current = iso_200['rmag'].values
color_current = g_current - r_current

# Estimate reddening from main sequence color at fixed mass
# For M = 0.5 M_sun (solar-type), expect (g-r)_0 ~ 0.45 mag
ms_mask = (iso_200['label'] == 0) & (iso_200['Mini'] > 0.45) & (iso_200['Mini'] < 0.55)
if np.any(ms_mask):
    color_obs = np.mean(color_current[ms_mask])
    color_intrinsic = 0.45  # Unreddened solar-type star
    E_gr_current = color_obs - color_intrinsic
    print(f"\n--- Estimated Current Reddening ---")
    print(f"Observed (g-r) for 0.5 M_sun MS star: {color_obs:.3f} mag")
    print(f"Expected intrinsic (g-r): {color_intrinsic:.3f} mag")
    print(f"Inferred E(g-r): {E_gr_current:.3f} mag")

    # For Cardelli+ 1989 with R_V = 3.1:
    # E(g-r) ≈ 0.98 × E(B-V)
    E_BV_current = E_gr_current / 0.98
    print(f"Inferred E(B-V): {E_BV_current:.3f}")
else:
    print("Warning: Could not find solar-type MS stars!")
    E_BV_current = 0.3  # Guess from visual inspection

# Cardelli+ 1989 extinction law (R_V = 3.1)
def extinction_coeffs(R_V=3.1):
    """
    Return A_λ/E(B-V) for SDSS filters.
    From Yuan+ 2013, Table 3 (Cardelli+ 1989 law)
    """
    # A_λ / E(B-V) for SDSS ugriz
    return {
        'u': 4.239,
        'g': 3.303,
        'r': 2.285,
        'i': 1.698,
        'z': 1.263
    }

coeffs = extinction_coeffs()
A_g_per_EBV = coeffs['g']
A_r_per_EBV = coeffs['r']

print(f"\nExtinction coefficients (Cardelli+ 1989, R_V=3.1):")
print(f"  A_g / E(B-V) = {A_g_per_EBV:.3f}")
print(f"  A_r / E(B-V) = {A_r_per_EBV:.3f}")

# Step 1: Remove current reddening to get intrinsic magnitudes at 10 pc
A_g_current = A_g_per_EBV * E_BV_current
A_r_current = A_r_per_EBV * E_BV_current

g_intrinsic_10pc = g_current - A_g_current
r_intrinsic_10pc = r_current - A_r_current
color_intrinsic = g_intrinsic_10pc - r_intrinsic_10pc

print(f"\n--- Step 1: Remove Current Reddening ---")
print(f"Removing A_g = {A_g_current:.3f} mag, A_r = {A_r_current:.3f} mag")
print(f"Intrinsic (g-r) at 10 pc: {color_intrinsic.min():.3f} - {color_intrinsic.max():.3f} mag")
print(f"Intrinsic g-band at 10 pc: {g_intrinsic_10pc.min():.2f} - {g_intrinsic_10pc.max():.2f} mag")

# Step 2: Apply distance modulus to M34 (470 pc)
distance_M34_pc = 470
mu_M34 = 5 * np.log10(distance_M34_pc / 10)  # From 10 pc to 470 pc

g_intrinsic_M34 = g_intrinsic_10pc + mu_M34
r_intrinsic_M34 = r_intrinsic_10pc + mu_M34

print(f"\n--- Step 2: Apply Distance Modulus ---")
print(f"Distance: 10 pc → {distance_M34_pc} pc")
print(f"Distance modulus: μ = {mu_M34:.3f} mag")
print(f"Intrinsic g-band at {distance_M34_pc} pc: {g_intrinsic_M34.min():.2f} - {g_intrinsic_M34.max():.2f} mag")

# Step 3: Apply M34 reddening
E_BV_M34 = 0.08
A_g_M34 = A_g_per_EBV * E_BV_M34
A_r_M34 = A_r_per_EBV * E_BV_M34

g_M34 = g_intrinsic_M34 + A_g_M34
r_M34 = r_intrinsic_M34 + A_r_M34
color_M34 = g_M34 - r_M34

print(f"\n--- Step 3: Apply M34 Extinction ---")
print(f"E(B-V) = {E_BV_M34:.3f}")
print(f"Adding A_g = {A_g_M34:.3f} mag, A_r = {A_r_M34:.3f} mag")
print(f"Final (g-r) range: {color_M34.min():.3f} - {color_M34.max():.3f} mag")
print(f"Final g-band range: {g_M34.min():.2f} - {g_M34.max():.2f} mag")

# Save corrected isochrone
iso_corrected = iso_200.copy()
iso_corrected['gmag'] = g_M34
iso_corrected['rmag'] = r_M34

output_file = 'isochrone_m34_corrected.dat'
iso_corrected.to_csv(output_file, sep=' ', index=False)
print(f"\n✓ Saved corrected isochrone to {output_file}")

# Visualization
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Panel 1: Original isochrone (10 pc, E(B-V)=0.93)
ax = axes[0, 0]
ax.plot(color_current, g_current, 'r-', lw=2, label=f'Original\n10 pc, E(B-V)≈{E_BV_current:.2f}')
ax.set_xlabel('g - r [mag]')
ax.set_ylabel('g [mag]')
ax.invert_yaxis()
ax.set_title(f'Step 0: Original PARSEC')
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_xlim(-0.5, 2.0)
ax.set_ylim(20, 0)

# Panel 2: Intrinsic at 10 pc (dereddened)
ax = axes[0, 1]
ax.plot(color_intrinsic, g_intrinsic_10pc, 'b-', lw=2, label='Intrinsic\n10 pc, E(B-V)=0')
ax.set_xlabel('g - r [mag]')
ax.set_ylabel('g [mag]')
ax.invert_yaxis()
ax.set_title('Step 1: Remove Reddening')
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_xlim(-0.5, 2.0)
ax.set_ylim(20, 0)

# Panel 3: Intrinsic at 470 pc (no reddening)
ax = axes[1, 0]
ax.plot(color_intrinsic, g_intrinsic_M34, 'purple', lw=2, label=f'Intrinsic\n{distance_M34_pc} pc, E(B-V)=0')
ax.set_xlabel('g - r [mag]')
ax.set_ylabel('g [mag]')
ax.invert_yaxis()
ax.set_title('Step 2: Apply Distance Modulus')
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_xlim(-0.5, 2.0)
ax.set_ylim(20, 0)

# Panel 4: Final M34 isochrone (470 pc, E(B-V)=0.08)
ax = axes[1, 1]
ax.plot(color_M34, g_M34, 'g-', lw=3, label=f'Final M34\n{distance_M34_pc} pc, E(B-V)={E_BV_M34}')
ax.set_xlabel('g - r [mag]')
ax.set_ylabel('g [mag]')
ax.invert_yaxis()
ax.set_title(f'Step 3: Apply M34 Extinction')
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_xlim(-0.5, 2.0)
ax.set_ylim(20, 0)

plt.tight_layout()
plt.savefig('isochrone_reddening_correction.png', dpi=150, bbox_inches='tight')
print(f"✓ Saved diagnostic plot to isochrone_reddening_correction.png")

print("\n" + "=" * 70)
print("CORRECTION COMPLETE!")
print("=" * 70)
print(f"\nCorrections Applied:")
print(f"  1. Reddening: E(B-V) = {E_BV_current:.2f} → {E_BV_M34:.2f}")
print(f"     - Removed A_g = {A_g_current:.2f} mag, A_r = {A_r_current:.2f} mag")
print(f"     - Added   A_g = {A_g_M34:.2f} mag, A_r = {A_r_M34:.2f} mag")
print(f"  2. Distance: 10 pc → {distance_M34_pc} pc")
print(f"     - Added distance modulus μ = {mu_M34:.2f} mag")
print(f"\nFinal Isochrone Properties:")
print(f"  Color range: (g-r) = {color_M34.min():.2f} to {color_M34.max():.2f} mag")
print(f"  Magnitude range: g = {g_M34.min():.2f} to {g_M34.max():.2f} mag")
print(f"  Expected M34 turnoff: g ~ 12-13 mag")
print(f"\nNext: Update membership.ipynb to use 'isochrone_m34_corrected.dat'")
print("=" * 70)
