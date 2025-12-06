#!/usr/bin/env python3
"""
Fix M2 spatial membership function

M2 is a globular cluster with ~150,000 members. Your image captures the
central ~20 arcmin, which is nearly all cluster members (not field stars!).

The background fraction needs to be estimated from the OUTER regions,
not set arbitrarily to 0.2.
"""

import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

print("=" * 70)
print("M2 SPATIAL MEMBERSHIP - PROPER BACKGROUND ESTIMATION")
print("=" * 70)

# Load catalog
catalog = pd.read_csv('m2_photometry_calibrated.csv')

# M2 parameters
ra_center = 323.3625
dec_center = -0.8231
r_core = 0.5  # arcmin

# Compute radial distances
coord_center = SkyCoord(ra=ra_center*u.deg, dec=dec_center*u.deg, frame='icrs')
coords_stars = SkyCoord(ra=catalog['ra'].values*u.deg, dec=catalog['dec'].values*u.deg, frame='icrs')
r_arcmin = coord_center.separation(coords_stars).arcmin

print(f"\nM2 photometry:")
print(f"  Total stars: {len(catalog)}")
print(f"  Radial range: {r_arcmin.min():.2f} to {r_arcmin.max():.2f} arcmin")

# Estimate background from OUTER region (r > 10 arcmin)
# This is beyond the tidal radius, so mostly field stars
outer_region = r_arcmin > 10.0  # arcmin

# Compute surface density in outer region
r_outer = r_arcmin[outer_region]
r_bins_outer = np.array([10.0, 15.0, 20.0, 30.0])
densities_outer = []

print(f"\n--- Background Estimation from Outer Region ---")
for i in range(len(r_bins_outer)-1):
    mask = (r_outer >= r_bins_outer[i]) & (r_outer < r_bins_outer[i+1])
    n = np.sum(mask)
    area = np.pi * (r_bins_outer[i+1]**2 - r_bins_outer[i]**2)
    density = n / area
    densities_outer.append(density)
    print(f"  {r_bins_outer[i]:.1f}-{r_bins_outer[i+1]:.1f} arcmin: {n} stars, σ = {density:.3f} stars/arcmin²")

# Background density (median of outer bins)
sigma_bg = np.median(densities_outer)
print(f"\nEstimated field background: σ_bg = {sigma_bg:.3f} stars/arcmin²")

# Compute center density
inner_region = r_arcmin < r_core
n_inner = np.sum(inner_region)
area_inner = np.pi * r_core**2
sigma_center = n_inner / area_inner

print(f"Center density (r < {r_core} arcmin): σ_0 = {sigma_center:.1f} stars/arcmin²")
print(f"Contrast: σ_0 / σ_bg = {sigma_center / sigma_bg:.0f}×")

# Background fraction relative to peak
bg_fraction = sigma_bg / sigma_center

print(f"\n--- New Background Fraction ---")
print(f"bg_fraction = σ_bg / σ_0 = {bg_fraction:.4f}")
print(f"  (This replaces the old value of 0.2)")
print(f"  → At center: P_spatial = {1.0 / (1.0 + bg_fraction):.3f}")
print(f"  → At r_core: {1.0 / (2**(5/2) + bg_fraction):.3f}")

# Compute spatial membership with proper background
def radial_membership_proper(r, r_core, bg_fraction):
    """
    Spatial membership using Plummer profile with properly estimated background.
    """
    rho_cluster = (1 + (r/r_core)**2)**(-2.5)
    rho_bg = bg_fraction
    prob = rho_cluster / (rho_cluster + rho_bg)
    return prob

P_spatial_new = radial_membership_proper(r_arcmin, r_core, bg_fraction)

print(f"\n--- Spatial Membership Statistics ---")
print(f"P_spatial range: {P_spatial_new.min():.4f} to {P_spatial_new.max():.4f}")
print(f"Median P_spatial: {np.median(P_spatial_new):.4f}")
print(f"Stars with P_spatial > 0.5: {np.sum(P_spatial_new > 0.5)} ({100*np.sum(P_spatial_new > 0.5)/len(P_spatial_new):.1f}%)")
print(f"Stars with P_spatial > 0.9: {np.sum(P_spatial_new > 0.9)} ({100*np.sum(P_spatial_new > 0.9)/len(P_spatial_new):.1f}%)")
print(f"Weighted spatial members: {np.sum(P_spatial_new):.1f} stars")

print("\n" + "=" * 70)
print("RECOMMENDATION FOR membership.ipynb")
print("=" * 70)
print(f"\nChange the background fraction calculation to:")
print(f"")
print(f"if cluster_name == 'M2':")
print(f"    # Estimate from outer region (r > 3 × r_core)")
print(f"    outer = r_arcmin > 3 * r_core")
print(f"    if np.sum(outer) > 10:")
print(f"        # Median density in outer region")
print(f"        r_bins_bg = np.array([3*r_core, 5*r_core, 10*r_core, 30*r_core])")
print(f"        densities = []")
print(f"        for i in range(len(r_bins_bg)-1):")
print(f"            mask = (r_arcmin >= r_bins_bg[i]) & (r_arcmin < r_bins_bg[i+1])")
print(f"            n = np.sum(mask)")
print(f"            area = np.pi * (r_bins_bg[i+1]**2 - r_bins_bg[i]**2)")
print(f"            densities.append(n / area)")
print(f"        sigma_bg = np.median(densities)")
print(f"        ")
print(f"        # Center density")
print(f"        inner = r_arcmin < r_core")
print(f"        sigma_center = np.sum(inner) / (np.pi * r_core**2)")
print(f"        ")
print(f"        bg_frac = sigma_bg / sigma_center")
print(f"    else:")
print(f"        bg_frac = 0.01  # Default for GC")
print(f"else:  # M34")
print(f"    bg_frac = 0.3  # Open cluster, more contamination")
print(f"")
print("=" * 70)
print(f"\nWith this change:")
print(f"  - M2 will have bg_fraction = {bg_fraction:.4f} (data-driven!)")
print(f"  - Spatial membership will be much higher at all radii")
print(f"  - Combined membership won't be crushed by low P_spatial")
print("=" * 70)
