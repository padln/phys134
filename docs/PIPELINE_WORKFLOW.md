# M34 Analysis Pipeline Workflow

## Architecture: .py modules + .ipynb visualization

**Philosophy**: Keep production code in `.py` files, use notebooks for debugging, visualization, and exploration.

---

## File Structure

```
PHYS 134L/
├── Data/
│   ├── tfn0m436-sq33-20251018-0206-e91.fits  (g-band)
│   └── tfn0m419-sq32-20251018-0154-e91.fits  (r-band)
│
├── pipeline_utils.py              ← Reusable photometry functions
├── data_reduction_simple.ipynb    ← Step 1: FITS → photometric catalog
├── membership.ipynb               ← Step 2: Catalog → cluster members
├── completeness.ipynb             ← Artificial star tests (parallel to Step 1)
├── density_profile.ipynb          ← Step 3: Members → radial profile
└── fisher_analysis.ipynb          ← Step 4: Observation optimization
```

---

## Complete Workflow

### Step 1: Data Reduction (`data_reduction_simple.ipynb`)

**Input**: Raw FITS files (g and r bands)

**Process**:
1. Load FITS images
2. Estimate and subtract 2D background
3. Detect sources (DAOStarFinder)
4. Aperture photometry with uncertainties
5. Cross-match g and r catalogs
6. Save catalog

**Output**: `m34_photometry_instrumental.csv`
- Columns: `ra`, `dec`, `x`, `y`, `g_inst`, `g_err`, `r_inst`, `r_err`, `g_minus_r`
- ~5000-8000 sources

**Uses**: Functions from `pipeline_utils.py`

---

### Step 2: Membership Determination (`membership.ipynb`)

**Input**: `m34_photometry_instrumental.csv` + Gaia data

**Process**:
1. Load photometric catalog
2. Query Gaia DR3 for proper motions
3. Cross-match catalog with Gaia
4. Load M34 isochrone (PARSEC/MIST)
5. **CMD filtering**: Distance to isochrone → P_CMD
6. **Proper motion filtering**: Cluster PM vs field → P_PM
7. **Spatial filtering**: Radial profile vs uniform → P_spatial
8. **Bayesian combination**: P_member = f(P_CMD, P_PM, P_spatial)
9. Filter catalog: keep stars with P_member > threshold

**Output**: `m34_members.csv`
- Same columns as input + `P_CMD`, `P_PM`, `P_spatial`, `P_member`
- ~200-400 likely members

**Key point**: This uses the **Bayesian membership framework** from your paper!

---

### Step 3: Photometric Calibration (in `membership.ipynb`)

**Process**:
1. Use Gaia photometry (G, BP, RP) for matched stars
2. Transform Gaia → SDSS g, r (color equations)
3. Calculate zeropoint: ZP = mag_Gaia - mag_inst
4. Apply calibration: g = g_inst + ZP_g, r = r_inst + ZP_r

**Output**: Add calibrated columns `g`, `r` to catalog

---

### Step 4: Completeness Correction (`completeness.ipynb`)

**Parallel to Step 1** - Run artificial star tests:

**Process**:
1. Add fake stars to FITS images at various magnitudes
2. Re-run detection + photometry
3. Measure recovery fraction C(m) vs magnitude
4. Fit completeness function (error function model)
5. Use Richardson-Lucy deconvolution for luminosity function

**Output**: Completeness function C(m) for both bands

**This feeds into Step 5** for density profile corrections

---

### Step 5: Density Profile Construction (`density_profile.ipynb`)

**Input**:
- `m34_members.csv` (cluster members only)
- Completeness function C(m)

**Process**:
1. Calculate distance from cluster center for each star
2. Bin stars in radial annuli
3. Apply completeness corrections: N_true(r) = N_obs(r) / C(m)
4. Calculate surface density: Σ(R) = N(R) / Area
5. Fit Plummer profile: Σ(R) = (M a²) / [π(R² + a²)²]
6. Extract parameters: core radius `a`, total mass `M`

**Output**:
- Radial profile plot
- Best-fit Plummer parameters
- Comparison with M2 (when you get M2 data)

---

### Step 6: Fisher Information Analysis (`fisher_analysis.ipynb`)

**Input**: Current photometry + completeness

**Process**:
1. Calculate Fisher information from existing data:
   - I = Σ (1/σ_m²) for all detected stars
   - Compute for Plummer parameters (a, M, μ)

2. Simulate additional observations:
   - Assume 2×, 3×, 5× exposure time
   - Predict new S/N and completeness
   - Calculate I_new

3. Diminishing returns analysis:
   - Plot I vs exposure time
   - Determine optimal observation strategy

**Output**:
- **Recommendation: Take more data? Which band? How long?**
- Quantitative answer to "is it worth it?"

---

## Data Flow Diagram

```
FITS files (g, r)
     ↓
[data_reduction_simple.ipynb]
     ↓
m34_photometry_instrumental.csv (~6000 sources)
     ↓
[membership.ipynb] + Gaia query + isochrone
     ↓
     ├→ P_CMD (distance to isochrone)
     ├→ P_PM (proper motion clustering)
     ├→ P_spatial (radial profile)
     └→ P_member = Bayesian combination
     ↓
m34_members.csv (~300 members, P_member > 0.5)
     ↓
[density_profile.ipynb] + completeness C(m)
     ↓
Radial profile + Plummer fit
     ↓
[fisher_analysis.ipynb]
     ↓
Observation recommendations
```

---

## Integration with Existing Code

### Your existing `membership.ipynb` already has:
✓ CMD isochrone distance calculation
✓ Proper motion clustering (iterative sigma clipping)
✓ Spatial distribution analysis
✓ Bayesian combination framework

**Just needs**: Real photometric catalog as input (instead of mock data)

### Your existing `completeness.ipynb` already has:
✓ Completeness function models (erf, tanh, Fermi-Dirac)
✓ Maximum likelihood fitting
✓ Richardson-Lucy deconvolution
✓ MCMC uncertainty estimation

**Just needs**: Real artificial star test data

---

## Running the Complete Analysis

### Phase 1: Get Photometry
```bash
# In Jupyter
1. Open data_reduction_simple.ipynb
2. Run all cells
3. Check output: m34_photometry_instrumental.csv
```

### Phase 2: Identify Members
```bash
4. Open membership.ipynb
5. Load m34_photometry_instrumental.csv
6. Query Gaia DR3
7. Apply Bayesian membership
8. Check output: m34_members.csv
```

### Phase 3: Build Profile
```bash
9. Run artificial star tests (completeness.ipynb)
10. Open density_profile.ipynb (needs to be created)
11. Load m34_members.csv
12. Construct corrected radial profile
13. Fit Plummer model
```

### Phase 4: Optimize Observations
```bash
14. Open fisher_analysis.ipynb
15. Calculate current Fisher information
16. Simulate additional exposures
17. Recommendation: more data needed? Yes/No
```

---

## Next Immediate Steps

1. **Test `data_reduction_simple.ipynb`** on your M34 data
   - Should take ~2 minutes to run
   - Check: Do you get ~6000 sources?
   - Check: Does CMD look reasonable?

2. **Modify `membership.ipynb`** to load CSV instead of generating mock data
   - Replace mock data generation with: `catalog = Table.read('m34_photometry_instrumental.csv')`
   - Add Gaia query for proper motions
   - Rest of code should work as-is!

3. **Create `density_profile.ipynb`** (I can help with this)
   - Read m34_members.csv
   - Bin by radius
   - Apply completeness corrections
   - Fit Plummer profile

4. **Run artificial star tests** with real FITS data
   - Inject fake stars into your images
   - Measure recovery fractions
   - Feed into completeness.ipynb

---

## Advantages of This Architecture

### Modular .py files:
✓ Reusable across projects
✓ Easy to unit test
✓ Version control friendly
✓ Can be imported anywhere

### Notebooks for visualization:
✓ See intermediate results
✓ Debug issues quickly
✓ Make plots for paper
✓ Explore parameter choices

### Clear data flow:
✓ Each step has defined input/output
✓ Can re-run any step independently
✓ Easy to share catalogs

---

## Questions?

- **Q: Do I need to calibrate photometry before membership?**
  A: No! Membership determination works fine with instrumental magnitudes. Calibrate after filtering to members only.

- **Q: Can I skip completeness corrections?**
  A: For membership? Yes. For density profile? No - completeness is critical for accurate mass/density.

- **Q: What if source catalogs don't match between g and r?**
  A: Cross-matching handles this. Unmatched sources are discarded (field star contamination).

- **Q: How do I know if my membership probabilities are good?**
  A: Check distributions: members should have P > 0.8, field stars P < 0.2, with clean separation.

---

Let me know if you want me to:
1. Create `density_profile.ipynb`
2. Create `fisher_analysis.ipynb`
3. Modify your existing `membership.ipynb` to work with this pipeline
4. Help debug the data reduction

**Test `data_reduction_simple.ipynb` first!**
