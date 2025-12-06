# M2 vs M34 Analysis Workflow

Complete guide for processing both clusters to compare density profiles

---

## 1. Isochrone Setup

### M34 (Open Cluster) - ALREADY DONE ✓
- **File:** `isochrone_m34_corrected.dat`
- Age: 200 Myr
- [M/H]: -0.1
- Distance: 470 pc
- E(B-V): 0.08

### M2 (Globular Cluster) - **NEEDS DOWNLOAD**

**Download from:** http://stev.oapd.inaf.it/cgi-bin/cmd

**Settings:**
```
Evolutionary tracks: PARSEC version 1.2S + COLIBRI
Photometric system: SDSS ugriz (Vega)
Circumstellar dust: OFF
Interstellar extinction: Av = 0 (apply manually)

Single isochrone:
  - log(age/yr) = 10.11  (13 Gyr)
  - [M/H] = -1.6  (or Z = 0.0001)
  - Initial mass range: 0.1 - 10 Msun

Output: Long format
```

**Save as:** `isochrone_m2_13gyr_mh-1.6.dat`

**Then apply corrections:**
```python
# M2 parameters
distance_M2 = 11500  # pc
mu_M2 = 5 * np.log10(distance_M2 / 10)  # 15.30 mag
E_BV_M2 = 0.06

# Apply distance modulus and reddening
g_M2 = g_parsec + mu_M2 + 3.303 * E_BV_M2
r_M2 = r_parsec + mu_M2 + 2.285 * E_BV_M2
```

---

## 2. Data Reduction Pipeline

### Status:
- ✅ **M34:** `data_reduction_simple.ipynb` (DONE)
- ✅ **M2:** `data_reduction_m2.ipynb` (CREATED - ready to run)

### Outputs:
- `m34_photometry_calibrated.csv`
- `m2_photometry_calibrated.csv` (after running M2 notebook)

---

## 3. Membership Determination

The `membership.ipynb` needs cluster-specific parameters. **Two options:**

### Option A: Parameterized Single Notebook (RECOMMENDED)

Add this cell at the top of `membership.ipynb`:

```python
# ============================================================================
# CLUSTER SELECTION - Change this to switch between M34 and M2
# ============================================================================

CLUSTER = "M34"  # Options: "M34" or "M2"

if CLUSTER == "M34":
    # M34 (Open Cluster) parameters
    cluster_name = "M34"
    catalog_file = "m34_photometry_calibrated.csv"
    isochrone_file = "isochrone_m34_corrected.dat"
    age_gyr = 0.2  # 200 Myr
    ra_center = 40.5092  # deg
    dec_center = 42.7619  # deg
    r_core = 5.0  # arcmin
    search_radius = 30.0 / 60.0  # deg

elif CLUSTER == "M2":
    # M2 (Globular Cluster) parameters
    cluster_name = "M2"
    catalog_file = "m2_photometry_calibrated.csv"
    isochrone_file = "isochrone_m2_13gyr_mh-1.6_corrected.dat"
    age_gyr = 13.0  # 13 Gyr
    ra_center = 323.3625  # deg (21h 33m 27s)
    dec_center = -0.8231  # deg (-0° 49' 23")
    r_core = 0.5  # arcmin (much smaller for globular!)
    search_radius = 30.0 / 60.0  # deg

print(f"Selected cluster: {cluster_name}")
print(f"  Catalog: {catalog_file}")
print(f"  Isochrone: {isochrone_file}")
print(f"  Center: RA={ra_center:.4f}°, Dec={dec_center:.4f}°")
```

Then replace all hardcoded values with these variables throughout the notebook.

### Option B: Two Separate Notebooks

```bash
cp membership.ipynb membership_m34.ipynb
cp membership.ipynb membership_m2.ipynb
```

Edit each with cluster-specific parameters.

---

## 4. Completeness Analysis

The `completeness.ipynb` notebook performs **artificial star tests** to determine detection completeness as a function of magnitude. This is **data-dependent**, so you need to run it separately for each cluster.

### Changes needed (minimal):

1. **Load different data:**
   ```python
   # For M34
   data_g = fits.open("./Data/tfn0m436-sq33-20251018-0206-e91.fits")[0].data

   # For M2
   data_g = fits.open("./Data/lco_data-20251126-4/lsc0m476-sq34-20251119-0090-e91.fits.fz")[1].data
   ```

2. **Change output filename:**
   ```python
   # Save results
   completeness_data.to_csv("m34_completeness.csv")  # or m2_completeness.csv
   ```

**Key difference:** M2 will have **much worse** completeness at faint magnitudes due to crowding (globular clusters are very dense!).

---

## 5. Mass Estimation

The `mass_estimation.ipynb` converts photometry to stellar masses using mass-luminosity relations and isochrones.

### Changes needed:

1. **Load appropriate isochrone:**
   ```python
   iso = load_parsec_isochrone(isochrone_file, age_gyr=age_gyr)
   ```

2. **Load member catalog (from membership analysis):**
   ```python
   members = pd.read_csv("m34_members.csv")  # or m2_members.csv
   ```

3. **Use cluster-specific parameters:**
   - M34: Solar metallicity, 200 Myr → main-sequence turnoff at ~2 M☉
   - M2: [Fe/H]=-1.6, 13 Gyr → turnoff at ~0.8 M☉

---

## 6. Density Profile Construction

This is the **final goal**! Create a new notebook: `density_profiles.ipynb`

### Workflow:

1. **Load member catalogs with probabilities:**
   ```python
   m34_members = pd.read_csv("m34_members.csv")  # Has P_member column
   m2_members = pd.read_csv("m2_members.csv")
   ```

2. **Compute radial distances:**
   ```python
   # Already computed in membership analysis, but recalculate if needed
   sep = coord_center.separation(coords_stars)
   r = sep.to(u.arcmin).value  # or u.pc after distance correction
   ```

3. **Apply completeness corrections:**
   ```python
   # Weight by completeness at each magnitude
   completeness = interp_completeness(members['g'])
   weights = members['P_member'] / completeness
   ```

4. **Bin by radius and compute density:**
   ```python
   r_bins = np.logspace(log10(r_min), log10(r_max), n_bins)

   for i in range(len(r_bins)-1):
       mask = (r >= r_bins[i]) & (r < r_bins[i+1])
       area = np.pi * (r_bins[i+1]**2 - r_bins[i]**2)  # arcmin²

       # Surface density (stars/arcmin²)
       sigma[i] = np.sum(weights[mask]) / area
       sigma_err[i] = np.sqrt(np.sum(weights[mask]**2)) / area
   ```

5. **Fit Plummer profile:**
   ```python
   def plummer_surface_density(R, M, a, sigma_bg):
       """Projected Plummer profile"""
       return (M * a**2) / (np.pi * (R**2 + a**2)**2) + sigma_bg

   # Fit
   popt, pcov = curve_fit(plummer_surface_density, r_centers, sigma,
                          sigma=sigma_err, p0=[1000, 5, 1])
   M_fit, a_fit, sigma_bg_fit = popt
   ```

6. **Plot side-by-side comparison:**
   ```python
   fig, axes = plt.subplots(1, 2, figsize=(14, 6))

   # M34
   ax = axes[0]
   ax.errorbar(r_m34, sigma_m34, yerr=sigma_err_m34, fmt='o', label='M34 data')
   ax.plot(r_fit, plummer(r_fit, *popt_m34), 'r-', label=f'Plummer fit (a={a_m34:.1f} arcmin)')
   ax.set_xlabel('Radius [arcmin]')
   ax.set_ylabel('Surface density [stars/arcmin²]')
   ax.set_xscale('log')
   ax.set_yscale('log')
   ax.set_title('M34 (Open Cluster)')

   # M2
   ax = axes[1]
   ax.errorbar(r_m2, sigma_m2, yerr=sigma_err_m2, fmt='s', label='M2 data')
   ax.plot(r_fit, plummer(r_fit, *popt_m2), 'b-', label=f'Plummer fit (a={a_m2:.1f} arcmin)')
   ax.set_xlabel('Radius [arcmin]')
   ax.set_ylabel('Surface density [stars/arcmin²]')
   ax.set_xscale('log')
   ax.set_yscale('log')
   ax.set_title('M2 (Globular Cluster)')
   ```

---

## 7. Expected Results

### M34 (Open Cluster):
- **Core radius:** ~5 arcmin (~0.7 pc at 470 pc)
- **Central density:** ~10-50 stars/arcmin²
- **Profile:** Shallow, extended
- **Total members:** ~200-400
- **Tidal radius:** ~15-20 arcmin (cluster is dispersing)

### M2 (Globular Cluster):
- **Core radius:** ~0.5 arcmin (~1.6 pc at 11.5 kpc)
- **Central density:** ~500-2000 stars/arcmin² (much denser!)
- **Profile:** Steep, concentrated
- **Total members:** ~10,000-50,000 (in our FOV)
- **Half-light radius:** ~6 arcmin (~20 pc)

### Key Comparison:
- M2 is **~100× denser** at the center
- M2 has **~10× smaller** core radius (in physical units)
- M34 profile will show **signs of disruption** (tidal truncation)
- M2 profile will show **strong central concentration** (well-fit by King/Plummer)

---

## 8. Action Items

### Immediate (do now):
1. ✅ Created `data_reduction_m2.ipynb`
2. ⬜ **Download M2 isochrone** from CMD web interface
3. ⬜ **Run** `data_reduction_m2.ipynb` to get `m2_photometry_calibrated.csv`
4. ⬜ Apply distance modulus and reddening to M2 isochrone

### Next (after data reduction):
5. ⬜ Parameterize `membership.ipynb` for both clusters
6. ⬜ Run membership analysis for M34 and M2
7. ⬜ Save member catalogs with probabilities

### Final (analysis):
8. ⬜ Run `completeness.ipynb` for both clusters
9. ⬜ Create `density_profiles.ipynb`
10. ⬜ Fit Plummer models and compare!

---

## 9. Key Python Snippets

### Distance conversion (arcmin → pc):
```python
# M34
d_M34 = 470  # pc
r_pc_M34 = r_arcmin * (d_M34 / 206265) * 60  # arcmin → pc

# M2
d_M2 = 11500  # pc
r_pc_M2 = r_arcmin * (d_M2 / 206265) * 60  # arcmin → pc
```

### Apply M2 isochrone corrections:
```python
# Create fix_isochrone_m2.py (similar to fix_isochrone_simple.py)
import pandas as pd
import numpy as np

iso = pd.read_csv('isochrone_m2_13gyr_mh-1.6.dat',
                   delim_whitespace=True, comment='#')

# M2 parameters
mu_M2 = 5 * np.log10(11500 / 10)  # 15.30 mag
E_BV = 0.06
A_g = 3.303 * E_BV  # 0.198 mag
A_r = 2.285 * E_BV  # 0.137 mag

# Apply corrections (PARSEC isochrone is at 10 pc with no reddening)
iso['gmag_M2'] = iso['gmag'] + mu_M2 + A_g
iso['rmag_M2'] = iso['rmag'] + mu_M2 + A_r

# Save
iso.to_csv('isochrone_m2_13gyr_mh-1.6_corrected.dat',
           sep=' ', index=False)
```

---

## 10. Common Pitfalls

1. **Units confusion:** M2 is at 11.5 kpc, M34 at 0.47 kpc → 24× farther!
2. **Crowding:** M2 will have severe crowding → completeness drops faster
3. **Field contamination:** M2 is toward Galactic center → more field stars
4. **Isochrone shape:** M2 (old) has red giant branch, M34 (young) does not
5. **Core radius:** Don't compare arcmin directly → convert to pc first!

---

## 11. Final Plot Ideas

1. **Side-by-side density profiles** (log-log)
2. **Cumulative mass profiles** M(< r)
3. **CMD comparison** with fitted isochrones
4. **Membership probability maps** (spatial)
5. **Completeness comparison** (M2 will be worse)

---

**Good luck! This is the core of your PHYS 134L project. 🌟**
