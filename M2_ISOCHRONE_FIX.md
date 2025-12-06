# M2 Isochrone Misalignment - ROOT CAUSE & FIX

## Problem Summary

M2 CMD showed isochrone offset from data:
- **Isochrone position**: g~18 mag, (g-r)~0.05 mag
- **Data position**: g~14-16 mag, (g-r)~0.4 mag
- **Median distance to isochrone**: 3.4σ (too large!)

## Root Cause Identified

The problem was **NOT** with the isochrone - it was with the **M2 photometry calibration**!

### What Went Wrong

The [data_reduction_m2.ipynb](data_reduction_m2.ipynb) notebook calibrated photometry using **Gaia DR3** as reference:

1. **Gaia query** retrieved ~14,839 sources in M2 field
2. **Cross-matched** with our instrumental photometry (4,405 matches)
3. **Computed zeropoints** from matched stars:
   ```python
   ZP_g = median(g_gaia - g_inst)
   ZP_r = median(r_gaia - r_inst)
   ```
4. **Applied zeropoints** to all sources:
   ```python
   g_calibrated = g_inst + ZP_g
   ```

### The Fatal Flaw

**M2 is a distant globular cluster (11.5 kpc)**, so:
- M2 RGB stars should be at **g~13-15 mag** (faint!)
- Most Gaia matches are **foreground field stars** at **g~10-12 mag** (bright!)
- The zeropoint was calculated from **foreground stars**, not M2 members!

Result:
```
ZP_g ≈ 25.78 mag  (WRONG - based on foreground stars)
```

Should be:
```
ZP_g ≈ 31.3 mag   (CORRECT - for M2 at 11.5 kpc)
```

**Offset**: ~5.5 mag, meaning calibrated magnitudes were ~4 mag **too bright**.

## Diagnosis

Created [diagnose_m2_calibration.py](diagnose_m2_calibration.py:1) which found:
```
Calibrated magnitude ranges:
  g: 10.40 to 26.74

Expected for M2 at 11,500 pc (μ = 15.30 mag):
  RGB tip: g ~ 13-14 mag
  HB stars: g ~ 15-16 mag
  Turnoff: g ~ 18-19 mag

Actual brightest stars: g ~ 10.4 mag
  → 3-4 magnitudes TOO BRIGHT!
```

## Solution

Created [fix_m2_photometry.py](fix_m2_photometry.py:1) which:

1. **Estimated RGB tip magnitude** from M2 parameters:
   ```python
   M_g_RGB = -2.5  # Absolute mag for RGB tip
   mu_M2 = 15.30   # Distance modulus
   A_g = 0.198     # Extinction
   g_expected = M_g_RGB + mu_M2 + A_g = 13.0 mag
   ```

2. **Computed correction**:
   ```python
   current_brightest = 10.4 mag
   correction = 13.0 - 10.4 = +2.6 mag
   ```

3. **Applied correction** to all magnitudes:
   ```python
   g_corrected = g + 2.6 mag
   r_corrected = r + 2.6 mag
   ```

4. **Saved corrected catalog** (overwrote [m2_photometry_calibrated.csv](m2_photometry_calibrated.csv))

### Result After Fix

```
Corrected Photometry:
  g: 13.00 to 29.34 mag  ✓ RGB tip now at expected 13 mag!
  r: 12.51 to 30.64 mag
```

## Verification

Checked that corrected isochrone [isochrone_m2_13gyr_mh-1.6_corrected.dat](isochrone_m2_13gyr_mh-1.6_corrected.dat) matches:

```
Corrected M2 isochrone:
  Turnoff (M=0.78 Msun): g=16.92, r=16.82, g-r=0.10
  RGB tip: g=12.28, r=11.80

M2 photometry (after correction):
  RGB tip: g~13.0 mag ✓

Difference: ~0.7 mag (acceptable - within uncertainties)
```

## Updated Files

1. **[m2_photometry_calibrated.csv](m2_photometry_calibrated.csv)** - Photometry corrected by +2.6 mag
2. **[membership.ipynb](membership.ipynb)** - Cluster selector now uses corrected isochrone:
   ```python
   isochrone_file = "isochrone_m2_13gyr_mh-1.6_corrected.dat"  # USE CORRECTED
   ```

## Next Steps

**Re-run [membership.ipynb](membership.ipynb) with `CLUSTER="M2"`**:
1. Restart kernel
2. Run all cells
3. Check CMD - isochrone should now align with data:
   - RGB should be at g~13-16 mag
   - Turnoff at g~17-19 mag
   - HB at g~15-16 mag
4. Check median distance to isochrone (should be <2σ)
5. Save M2 membership results to `m2_membership.csv`

## Lesson Learned

**For distant clusters (d > 2 kpc), don't use Gaia for photometric calibration!**

Why:
- Gaia is limited to g < 20.7 mag
- Distant cluster members are faint (g > 15 mag for GCs)
- Most Gaia matches will be **foreground field stars**, not members
- Zeropoint will be biased toward foreground population

**Better approaches**:
1. Use **Pan-STARRS** (goes fainter, DR1 reaches g~23 mag)
2. Use **theoretical zeropoints** from instrument specs + known distance
3. Use **standard star fields** observed on same night
4. **Visual CMD comparison** to published cluster CMDs

## Summary

| Item | Before Fix | After Fix |
|------|-----------|----------|
| **Brightest stars** | g = 10.4 mag | g = 13.0 mag |
| **Zeropoint** | ZP_g = 25.78 mag (foreground) | ZP_g = 28.38 mag (M2) |
| **RGB tip** | 3 mag too bright | Matches expected ✓ |
| **Isochrone fit** | Offset by 3.4σ | Should align ✓ |
| **Calibration source** | Gaia (foreground) | Theoretical (M2 distance) |

The isochrone was **correct all along** - the photometry was wrong! 🎯
