# IMPLEMENTATION CHECKLIST: Completing the Research Paper

**Goal:** Transform `article.tex` from framework → publication-ready manuscript  
**Timeline:** Sequential steps, each ~30-60 min  
**Outcome:** Submission-ready PDF with all figures, tables, and final results

---

## PHASE 1: VERIFY CURRENT CONTENT (15 min)

- [ ] Open `article.tex` in editor
- [ ] Compile with `pdflatex article.tex && bibtex article && pdflatex article.tex` (twice)
- [ ] Verify PDF generation without errors (ignore missing figure warnings for now)
- [ ] Review PDF page count (expect ~15-18 pages with 6 figures)
- [ ] Check that all sections render correctly

**Expected Issues & Fixes:**
- Missing `.png` files → Expected, will add in Phase 2
- Undefined references → Normal until figures added
- Citation warnings → Will resolve in Phase 3

---

## PHASE 2: GENERATE FIGURES FROM NOTEBOOKS (90 min)

### 2.1 CMD with Membership Probabilities (20 min)

**From:** `membership.ipynb` (wherever CMD plot cells are)  
**Produce:** `cmd_membership_m34.png` and `cmd_membership_m2.png`

**Steps:**
```python
# In membership.ipynb, add this cell after CMD plotting:

import matplotlib.pyplot as plt

# FIGURE: M34 CMD with membership colors
fig, ax = plt.subplots(figsize=(8, 10))
scatter = ax.scatter(
    catalog['g_minus_r'], catalog['g'],
    c=P_cmd,
    s=20,
    cmap='coolwarm',
    vmin=0,
    vmax=1,
    alpha=0.6,
    edgecolors='none'
)
# Overlay isochrone
ax.plot(color_iso_m34, g_iso_m34, 'k-', linewidth=2, label='PARSEC isochrone (200 Myr, [M/H]=0.0)')
ax.set_xlabel('$g\' - r\'$ (mag)', fontsize=12)
ax.set_ylabel('$g\'$ (mag)', fontsize=12)
ax.set_title('M34: Color-Magnitude Diagram with Membership', fontsize=14, fontweight='bold')
ax.invert_yaxis()
ax.grid(alpha=0.3)
cbar = plt.colorbar(scatter, ax=ax)
cbar.set_label('$P_{CMD}$ (Membership Probability)', fontsize=11)
ax.legend(loc='upper left', fontsize=11)
plt.tight_layout()
plt.savefig('../AASTeX_Template/cmd_membership_m34.png', dpi=300, bbox_inches='tight')
plt.close()

# Repeat for M2 with M2 isochrone data
# ... [similar code for M2]
plt.savefig('../AASTeX_Template/cmd_membership_m2.png', dpi=300, bbox_inches='tight')
```

**Verification:**
- [ ] File size ~500 KB-2 MB (acceptable for AASTeX)
- [ ] Isochrone clearly visible
- [ ] Color gradient shows membership variation
- [ ] Caption ready: "Color-magnitude diagrams for M34 and M2 with CMD membership probabilities..."

---

### 2.2 Spatial Distribution Maps (25 min)

**From:** `membership.ipynb` (spatial distribution plotting cells)  
**Produce:** `spatial_distribution_m34.png` and `spatial_distribution_m2.png`

**Steps:**
```python
# Two-panel plot showing both clusters
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

# M34
scatter1 = ax1.scatter(
    catalog_m34['ra'], catalog_m34['dec'],
    s=50 * P_combined_m34,  # Size proportional to membership
    c=P_combined_m34,
    cmap='viridis',
    vmin=0,
    vmax=1,
    alpha=0.6,
    edgecolors='k',
    linewidths=0.5
)
ax1.plot(ra_center_m34, dec_center_m34, 'r+', markersize=20, markeredgewidth=2)
ax1.set_xlabel('RA (deg)', fontsize=12)
ax1.set_ylabel('Dec (deg)', fontsize=12)
ax1.set_title('M34: Spatial Distribution', fontsize=14, fontweight='bold')
ax1.grid(alpha=0.3)

# M2 (similar)
scatter2 = ax2.scatter(
    catalog_m2['ra'], catalog_m2['dec'],
    s=50 * P_combined_m2,
    c=P_combined_m2,
    cmap='viridis',
    vmin=0,
    vmax=1,
    alpha=0.6,
    edgecolors='k',
    linewidths=0.5
)
ax2.plot(ra_center_m2, dec_center_m2, 'r+', markersize=20, markeredgewidth=2)
ax2.set_xlabel('RA (deg)', fontsize=12)
ax2.set_ylabel('Dec (deg)', fontsize=12)
ax2.set_title('M2: Spatial Distribution', fontsize=14, fontweight='bold')
ax2.grid(alpha=0.3)

cbar = plt.colorbar(scatter1, ax=ax2, label='$P_{combined}$')
plt.suptitle('Cluster Membership Probabilities by Sky Position', fontsize=15, fontweight='bold', y=1.00)
plt.tight_layout()
plt.savefig('../AASTeX_Template/spatial_distribution_both.png', dpi=300, bbox_inches='tight')
plt.close()
```

**Verification:**
- [ ] Both clusters visible with clear central concentration
- [ ] Color gradient visible (bright = high membership)
- [ ] Cluster centers marked with crosses
- [ ] File size acceptable

---

### 2.3 Radial Density Profiles (30 min)

**From:** `membership.ipynb` (radial profile cells)  
**Produce:** `radial_profiles_both.png`

**Steps:**
```python
# Already mostly done in membership.ipynb, just export high-quality version
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# M34 profile
ax1.errorbar(r_centers_m34, sigma_m34, yerr=sigma_err_m34,
             fmt='o', color='blue', markersize=8, capsize=4, capthick=1.5, 
             elinewidth=1.5, label='Weighted members')
ax1.axhline(bg_m34, color='gray', linestyle='--', alpha=0.7, label='Field background')
ax1.loglog()
ax1.set_xlabel('Radius (arcmin)', fontsize=12)
ax1.set_ylabel('Surface Density (stars/arcmin²)', fontsize=12)
ax1.set_title('M34: Radial Profile', fontsize=13, fontweight='bold')
ax1.grid(alpha=0.3, which='both')
ax1.legend()

# M2 profile (similar)
ax2.errorbar(r_centers_m2, sigma_m2, yerr=sigma_err_m2,
             fmt='s', color='red', markersize=8, capsize=4, capthick=1.5,
             elinewidth=1.5, label='Weighted members')
ax2.axhline(bg_m2, color='gray', linestyle='--', alpha=0.7, label='Field background')
ax2.loglog()
ax2.set_xlabel('Radius (arcmin)', fontsize=12)
ax2.set_ylabel('Surface Density (stars/arcmin²)', fontsize=12)
ax2.set_title('M2: Radial Profile', fontsize=13, fontweight='bold')
ax2.grid(alpha=0.3, which='both')
ax2.legend()

plt.suptitle('Radial Surface Density Profiles', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('../AASTeX_Template/radial_profiles_both.png', dpi=300, bbox_inches='tight')
plt.close()
```

**Verification:**
- [ ] Log-log scale shows power-law behavior clearly
- [ ] Error bars visible and reasonable
- [ ] Background levels clearly marked
- [ ] Power-law slopes estimable by eye

---

### 2.4 Background Estimation Panels (15 min)

**From:** Existing files or re-generate if needed  
**Status:** May already exist (`background_panel.png`, `background_panel_m34.png`)

**Check:**
- [ ] File `AASTeX_Template/background_panel.png` exists (3-panel M2)
- [ ] File `AASTeX_Template/background_panel_m34.png` exists (3-panel M34)
- [ ] Show: original image | background model | subtracted image

If missing, regenerate:
```python
from photutils.background import Background2D, MedianBackground
from photutils.detection import DAOStarFinder

# Load FITS image (for M34)
data = fits.open('Data/m34_g_image.fits')[0].data

# Estimate 2D background
bkg = Background2D(data, box_size=64, filter_size=3, bkg_estimator=MedianBackground())
bkg_image = bkg.background
data_sub = data - bkg_image

# Create 3-panel figure
fig = plt.figure(figsize=(15, 4))
ax1, ax2, ax3 = fig.subplots(1, 3)

im1 = ax1.imshow(data, cmap='gray', norm=PercentileNorm(vmin=1, vmax=99))
ax1.set_title('Original Image')
plt.colorbar(im1, ax=ax1)

im2 = ax2.imshow(bkg_image, cmap='gray')
ax2.set_title('2D Background Model')
plt.colorbar(im2, ax=ax2)

im3 = ax3.imshow(data_sub, cmap='gray', norm=PercentileNorm(vmin=1, vmax=99))
ax3.set_title('Background-Subtracted')
plt.colorbar(im3, ax=ax3)

plt.suptitle('Background Estimation: M34 g-band', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('../AASTeX_Template/background_panel_m34.png', dpi=200, bbox_inches='tight')
plt.close()
```

---

## PHASE 3: UPDATE ARTICLE.TEX WITH FINAL NUMBERS (30 min)

### 3.1 Verify Table Values Match Code

**Table 1 (Cluster Properties):**
- [ ] Distance values from literature
- [ ] Age values confirmed
- [ ] [Fe/H] values correct
- [ ] Concentration parameter $c$ recalculated if different

**Table 2 (Membership Statistics):**
- [ ] Re-run: `sum(P_combined > 0.5)` for each cluster
- [ ] Re-run: `sum(P_combined > 0.7)` for each cluster
- [ ] Re-run: `sum(P_combined > 0.9)` for each cluster
- [ ] Update percentages: `(count / total) * 100`
- [ ] Update effective counts: `sum(P_combined)`

**Script to extract values:**
```python
import pandas as pd
members = pd.read_csv('m34_membership.csv')

# M34
print(f"M34 Membership Statistics:")
print(f"Total sources: {len(members)}")
print(f"P > 0.1: {sum(members['P_combined'] > 0.1)} ({100*sum(members['P_combined'] > 0.1)/len(members):.1f}%)")
print(f"P > 0.3: {sum(members['P_combined'] > 0.3)} ({100*sum(members['P_combined'] > 0.3)/len(members):.1f}%)")
print(f"P > 0.5: {sum(members['P_combined'] > 0.5)} ({100*sum(members['P_combined'] > 0.5)/len(members):.1f}%)")
print(f"P > 0.7: {sum(members['P_combined'] > 0.7)} ({100*sum(members['P_combined'] > 0.7)/len(members):.1f}%)")
print(f"P > 0.9: {sum(members['P_combined'] > 0.9)} ({100*sum(members['P_combined'] > 0.9)/len(members):.1f}%)")
print(f"Sum of probabilities: {sum(members['P_combined']):.1f}")
print(f"Mean P (weighted): {sum(members['P_combined']) / len(members['P_combined']>0):.3f}")
```

**Table 3 (Plummer Fits):**
- [ ] $a$ values with uncertainties from fit output
- [ ] $M$ values from best-fit parameters
- [ ] $\chi^2_\nu$ from goodness-of-fit
- [ ] Literature values: double-check citations

### 3.2 Update Results Section Narrative

Search article.tex for placeholder text like `[to be updated]` and replace with actual values:
- [ ] M34 member count: 183 → actual count
- [ ] M2 member count: 2,585 → actual count
- [ ] M2 Plummer radius: 1.2 arcmin → actual fit
- [ ] M2 total mass: $1.2 \times 10^4$ → actual fit
- [ ] M34 Plummer radius: 6.5 arcmin → actual fit
- [ ] M34 total mass: $8.0 \times 10^2$ → actual fit

---

## PHASE 4: COMPILE AND VALIDATE (20 min)

### 4.1 Final Compilation
```bash
cd AASTeX_Template/

# Clean previous build
rm -f *.aux *.bbl *.blg *.log *.out

# Compile (3-pass to ensure all references resolved)
pdflatex -interaction=nonstopmode article.tex
bibtex article
pdflatex -interaction=nonstopmode article.tex
pdflatex -interaction=nonstopmode article.tex
```

**Expected output:** `article.pdf` without errors

### 4.2 PDF Quality Check
- [ ] All pages render (expect ~17-20 pages)
- [ ] All figures display clearly
- [ ] All table values readable
- [ ] Bibliography appears correctly at end
- [ ] No [??] references or ???
- [ ] Equations are numbered and formatted properly
- [ ] Font sizes are consistent

### 4.3 Specific Checks
- [ ] Page 1: Title, authors, abstract
- [ ] Page 2-3: Introduction with Table 1
- [ ] Page 4-7: Methodology sections with equations
- [ ] Page 8-9: Results with Table 2 and figures
- [ ] Page 9-10: Plummer fits (Table 3) and completeness assessment
- [ ] Page 11-13: Discussion with dynamical evolution details
- [ ] Page 14: Conclusions
- [ ] Page 15+: Appendices with figures

---

## PHASE 5: CITATION & REFERENCE AUDIT (15 min)

### 5.1 Check All Citations Exist in references.bib
```bash
# Get list of \cite commands
grep -o '\\cite[t]*{[^}]*}' article.tex | sort -u

# Verify each key exists in references.bib
```

**Required citations to verify:**
- [x] Plummer 1911
- [x] Harris 1996
- [x] Henry & McCarthy 2004
- [x] Kroupa 2001
- [x] PARSEC (Bressan 2012)
- [x] Gaia DR3
- [x] Richardson 1972, Lucy 1974
- [x] Photutils, AstroPy

### 5.2 Add Missing Keys
If any \cite command fails, add to references.bib:
```bibtex
@article{reference_key,
  author = {Author, A. and Name, B.},
  year = {2020},
  title = {Paper Title},
  journal = {Journal},
  volume = {XX},
  pages = {YY-ZZ}
}
```

---

## PHASE 6: SUPPLEMENTARY MATERIALS (Optional, 30 min)

### 6.1 Create Source Catalogs (CSV)
```python
# Export membership catalog for M34
members_m34 = pd.read_csv('m34_membership.csv')
members_m34[['ra', 'dec', 'g', 'r', 'P_cmd', 'P_spatial', 'P_combined']].to_csv(
    'AASTeX_Template/m34_members_catalog.csv',
    index=False,
    float_format='%.6f'
)
print(f"Saved {len(members_m34)} sources to m34_members_catalog.csv")

# Same for M2
```

### 6.2 Create Supplementary Table: Completeness Results
Once AST is run, create table:
```
Magnitude | N_added | N_recovered | C(m)    | σ_C
----------|---------|-------------|---------|--------
18.0      | 100     | 98          | 0.98    | 0.014
18.5      | 100     | 96          | 0.96    | 0.019
...
```

### 6.3 Create Appendix D (Optional): Full IMF Table
```
Mass (M☉) | N_observed | N_corrected | ξ(M)     | Error
----------|-----------|-------------|----------|--------
0.10      | 0         | 25          | 15.3     | ±2.1
0.15      | 2         | 45          | 28.5     | ±4.2
...
```

---

## FINAL CHECKLIST: SUBMISSION-READY

- [ ] PDF compiles without errors or warnings
- [ ] All figures display (6 required PNG files present)
- [ ] All tables have current data
- [ ] All citations resolved (no ???)
- [ ] Bibliography formatted correctly
- [ ] Word count ~8,500 (acceptable for AJ/ApJ)
- [ ] Equations numbered and referenced properly
- [ ] Author names and affiliations correct
- [ ] Abstract is self-contained and ~250 words
- [ ] Keywords present and UAT-aligned
- [ ] Acknowledgments mention LCO, Gaia
- [ ] Software cited (Photutils, AstroPy, emcee)

---

## SUBMISSION VENUES & FORMAT

### Recommended: Astronomical Journal (AJ)
- Format: AASTeX 7.0 (current template)
- Max length: ~25 pages (we have ~18, acceptable)
- Submission: arXiv first (recommended)

### Alternative: Astrophysical Journal (ApJ)
- Same AASTeX format
- Slightly stricter on methodology rigor
- Good fit for this work

### Process:
1. Create account on arXiv.org
2. Upload article.pdf, references.bib, and all PNG figures
3. Review preview, submit to astronomy/high-energy astrophysics
4. Wait for email confirmation
5. Once on arXiv, submit to AJ via the journal's editorial system
6. Include arXiv ID (arXiv:XXXX.XXXXX) in journal submission

---

## TIME ESTIMATE

| Phase | Task | Time |
|-------|------|------|
| 1 | Content verification | 15 min |
| 2 | Figure generation | 90 min |
| 3 | Update values/numbers | 30 min |
| 4 | Compile & validate | 20 min |
| 5 | Citation audit | 15 min |
| 6 | Supplementary materials | 30 min |
| **TOTAL** | | **200 min (3.3 hours)** |

---

## Questions & Troubleshooting

**Q: pdflatex fails with "file not found"**  
A: Ensure PNG files are in `AASTeX_Template/` directory, same as article.tex

**Q: Figure appears but is too small/large**  
A: Adjust `width=0.48\textwidth` in `\includegraphics{..}` command

**Q: Bibliography shows [?]**  
A: Run `bibtex article` and recompile 3 times

**Q: Table values don't match my analysis**  
A: Re-run Python scripts above to extract correct values, update by hand in article.tex

**Q: Can't find references.bib**  
A: Should be in `AASTeX_Template/` directory; if missing, retrieve from repo

---

This completes the implementation guide for finalizing the research paper!
