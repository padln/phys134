#!/usr/bin/env python3
"""
Apply distance modulus and reddening corrections to M2 isochrone

M2 (NGC 7089) parameters:
- Distance: 11,500 pc (distance modulus μ = 15.30 mag)
- E(B-V) = 0.06 mag
- Age: 13 Gyr
- [M/H] = -1.6

PARSEC isochrones are at 10 pc with no reddening, so we need to:
1. Add distance modulus μ = 5*log10(d/10 pc)
2. Add reddening: A_g = 3.303 * E(B-V), A_r = 2.285 * E(B-V)
"""

import pandas as pd
import numpy as np

# Input/output files
input_file = "isochrone_m2_13gyr_mh-1.6.dat"
output_file = "isochrone_m2_13gyr_mh-1.6_corrected.dat"

print("=" * 70)
print("M2 ISOCHRONE CORRECTION")
print("=" * 70)

# Load PARSEC isochrone
print(f"\nLoading isochrone from: {input_file}")

# Column names from PARSEC
column_names = [
    'Zini', 'MH', 'logAge', 'Mini', 'int_IMF', 'Mass', 'logL', 'logTe',
    'logg', 'label', 'McoreTP', 'C_O', 'period0', 'period1', 'period2',
    'period3', 'period4', 'pmode', 'Mloss', 'tau1m', 'X', 'Y', 'Xc',
    'Xn', 'Xo', 'Cexcess', 'Z', 'mbolmag', 'umag', 'gmag', 'rmag',
    'imag', 'zmag'
]

iso = pd.read_csv(
    input_file,
    delim_whitespace=True,
    comment='#',
    names=column_names,
    na_values=['-1.00', '-9.999']
)

# Clean up any header artifacts
for col in column_names:
    iso[col] = pd.to_numeric(iso[col], errors='coerce')
iso = iso.dropna(subset=['logAge', 'Mini']).reset_index(drop=True)

print(f"✓ Loaded {len(iso)} isochrone points")
print(f"  Age: {10**(iso['logAge'].iloc[0]-9):.1f} Gyr")
print(f"  [M/H]: {iso['MH'].iloc[0]:.2f}")
print(f"  Mass range: {iso['Mini'].min():.2f} - {iso['Mini'].max():.1f} M☉")

# M2 parameters
distance_M2 = 11500  # pc
mu_M2 = 5 * np.log10(distance_M2 / 10)  # distance modulus
E_BV_M2 = 0.06  # E(B-V)

# Reddening (Cardelli+ 1989, R_V = 3.1)
A_g = 3.303 * E_BV_M2
A_r = 2.285 * E_BV_M2

print(f"\n--- M2 Parameters ---")
print(f"Distance: {distance_M2:,} pc")
print(f"Distance modulus: μ = {mu_M2:.2f} mag")
print(f"Reddening: E(B-V) = {E_BV_M2:.3f} mag")
print(f"Extinction: A_g = {A_g:.3f} mag, A_r = {A_r:.3f} mag")

# Apply corrections
# PARSEC magnitudes are at 10 pc with no reddening
g_orig = iso['gmag'].values
r_orig = iso['rmag'].values

g_M2 = g_orig + mu_M2 + A_g
r_M2 = r_orig + mu_M2 + A_r

print(f"\n--- Magnitude Corrections ---")
print(f"g-band: {g_orig[0]:.2f} → {g_M2[0]:.2f} (Δ = {g_M2[0]-g_orig[0]:.2f} mag)")
print(f"r-band: {r_orig[0]:.2f} → {r_M2[0]:.2f} (Δ = {r_M2[0]-r_orig[0]:.2f} mag)")

# Create corrected isochrone
iso_corrected = iso.copy()
iso_corrected['gmag'] = g_M2
iso_corrected['rmag'] = r_M2
iso_corrected['imag'] = iso['imag'] + mu_M2 + 1.569 * E_BV_M2  # A_i
iso_corrected['zmag'] = iso['zmag'] + mu_M2 + 1.196 * E_BV_M2  # A_z
iso_corrected['umag'] = iso['umag'] + mu_M2 + 4.239 * E_BV_M2  # A_u

# Save
iso_corrected.to_csv(output_file, sep=' ', index=False, float_format='%.6f')
print(f"\n✓ Saved corrected isochrone to: {output_file}")

# Verification
print(f"\n--- Verification ---")
print(f"Select a solar-type star (M ~ 0.8 M☉, near turnoff):")
solar_idx = np.argmin(np.abs(iso['Mini'] - 0.8))
print(f"  Mass: {iso['Mini'].iloc[solar_idx]:.2f} M☉")
print(f"  Original: g={g_orig[solar_idx]:.2f}, r={r_orig[solar_idx]:.2f}, (g-r)={g_orig[solar_idx]-r_orig[solar_idx]:.2f}")
print(f"  Corrected: g={g_M2[solar_idx]:.2f}, r={r_M2[solar_idx]:.2f}, (g-r)={g_M2[solar_idx]-r_M2[solar_idx]:.2f}")
print(f"\n  Expected for M2 turnoff:")
print(f"    g ~ 18-19 mag (distance modulus + ~3 mag absolute)")
print(f"    (g-r) ~ 0.4-0.5 mag (metal-poor, old)")

# Check evolutionary phases
print(f"\n--- Evolutionary Phases ---")
phase_names = {
    0: 'Main Sequence',
    1: 'Subgiant',
    2: 'RGB base',
    3: 'RGB tip',
    4: 'HB',
    5: 'Early AGB',
    6: 'TP-AGB',
    7: 'post-AGB',
    8: 'WR',
    9: 'White Dwarf'
}

unique_labels = np.unique(iso['label'].values)
for label in unique_labels:
    n = np.sum(iso['label'] == label)
    phase = phase_names.get(label, f'Unknown ({label})')
    print(f"  {phase}: {n} points")

print("\n" + "=" * 70)
print("COMPLETE!")
print("=" * 70)
print(f"\nNext steps:")
print(f"  1. Use '{output_file}' in membership.ipynb")
print(f"  2. Load with load_parsec_isochrone('{output_file}', age_gyr=13.0)")
print(f"  3. M2 will have prominent RGB (evolved stars) - filter if needed")
print("=" * 70)
