# DETAILED IMPLEMENTATION CHECKLIST
## Stellar Cluster Density Profile Study - Gap Completion Guide

---

## NOTEBOOK-BY-NOTEBOOK ANALYSIS

### 1. signal_noise.ipynb - Status: 40% Complete

#### Current Functionality:
```python
✓ Basic S/N calculation (Equation 1)
✓ Multi-pixel aperture S/N (Equation 2)
✓ S/N vs. integration time plots
✓ Example calculations for M2 star
```

#### Missing Components:
```
❌ Exposure time calculator function
❌ Trade-off analysis: S/N vs. saturation vs. integration time
❌ Crowding penalty function
❌ Stacking strategy optimizer (n_exp, τ combinations)
❌ LCO telescope specification integration
❌ Atmospheric extinction effects (airmass dependence)
❌ Detector gain/readout noise for actual LCO cameras
❌ Filter transmission curves
❌ Observation planning wizard
```

#### Enhancement Tasks:
- [ ] Add function: `exposure_time_calculator(magnitude, target_snr, ...)`
- [ ] Add function: `stacking_optimization(mag_range, available_time, ...)`
- [ ] Query LCO API for actual telescope specs
- [ ] Plot S/N requirements vs. observing conditions
- [ ] Create interactive S/N calculator
- [ ] Validate against published LCO sensitivity curves

---

### 2. completeness.ipynb - Status: 80% Complete

#### Current Functionality:
```python
✓ Error function completeness model (Equation 12)
✓ Binomial log-likelihood (Equations 8-10)
✓ MLE fitting with Nelder-Mead optimization
✓ Hessian-based Fisher information
✓ MCMC sampling (emcee)
✓ Richardson-Lucy deconvolution (Equations 15-17)
✓ Corner plots and posterior visualization
```

#### Missing Components:
```
❌ Real artificial star test data integration
❌ Hyperbolic tangent model (Equation 13)
❌ Fermi-Dirac model (Equation 14)
❌ Model comparison/selection (AIC/BIC)
❌ Spatial completeness variation
❌ Magnitude-dependent systematics
❌ Color-dependent completeness
❌ Bootstrap resampling for systematic errors
❌ Deconvolution convergence diagnostics
❌ Comparison to simple C(m) = N_rec/N_add
```

#### Enhancement Tasks:
- [ ] Implement Equations 13-14 completeness models
- [ ] Add model comparison: `compare_completeness_models()`
- [ ] Create spatial completeness map function
- [ ] Implement `bootstrap_completeness_errors()`
- [ ] Add deconvolution convergence plots
- [ ] Write function: `apply_completeness_correction(catalog, C_func, ...)`
- [ ] Test on synthetic luminosity functions
- [ ] Add docstrings and parameter descriptions
- [ ] Create validation plots (residuals, predictions)

---

### 3. membership.ipynb - Status: 70% Complete

#### Current Functionality:
```python
✓ CMD membership: isochrone distance calculation
✓ Proper motion determination (sigma-clipping iteration)
✓ PM membership: Gaussian likelihood
✓ Spatial membership: radial profile weighting
✓ Combined Bayesian membership (Equation 25)
✓ Vector point diagram visualization
✓ Multi-panel membership visualization
```

#### Missing Components:
```
❌ Real Gaia DR3 data (uses mock)
❌ Real PARSEC isochrones (uses mock)
❌ Iterative membership refinement loop
❌ Field contamination background model
❌ Independence assumption validation
❌ Sensitivity to prior choice
❌ Alternative methods:
   ❌ Hierarchical Bayesian mixture model
   ❌ Machine learning classifier (Random Forest)
   ❌ Gaussian mixture model (GMM)
❌ Cross-validation of membership
❌ Comparison: PM-only vs. CMD-only vs. spatial-only
❌ Magnitude/color dependence of membership quality
```

#### Enhancement Tasks:
- [ ] Download actual M2 & M34 Gaia DR3 data: `gaia_data = fetch_gaia_cluster(...)`
- [ ] Download PARSEC isochrones for M2 & M34 parameters
- [ ] Implement iterative refinement loop
- [ ] Add: `hierarchical_bayesian_membership()` using PyMC3
- [ ] Add: `random_forest_membership()` for validation
- [ ] Implement cross-validation: `leave_one_out_cv()`
- [ ] Create sensitivity analysis for prior_member parameter
- [ ] Quantify field contamination
- [ ] Generate membership quality metrics
- [ ] Compare method combinations (test all 7 possible combinations)

---

### 4. mass_estimation.ipynb - Status: 50% Complete

#### Current Functionality:
```python
✓ Henry & McCarthy MLR for low-mass stars
✓ Mock isochrone generation (age/[Fe/H] dependent)
✓ Isochrone interpolator class
✓ Kroupa IMF implementation
✓ Binary fraction correction framework
✓ MCMC mass inference setup
✓ Plummer profile fitting
✓ Plummer surface density model
```

#### Missing Components:
```
❌ Real PARSEC/MIST isochrones loaded
❌ 2D isochrone interpolation (M,g,g-r)
❌ Extinct stars handling
❌ RGB/AGB star treatment
❌ Main sequence turnoff identification
❌ Realistic MCMC likelihood calibration
❌ Individual star mass posterior sampling
❌ Population-level IMF fitting
❌ Binary fraction inference
❌ Mass segregation analysis
❌ Extinction law application (Cardelli)
❌ Parallax/distance modulus handling
❌ Comparison of MLR vs. isochrone vs. spectroscopy
```

#### Enhancement Tasks:
- [ ] Write: `download_parsec_isochrone(age_gyr, feh)`
- [ ] Implement proper 2D interpolation: `Isochrone.mass_from_colors_mag()`
- [ ] Add extinction law: `cardelli_extinction(wavelength, Av)`
- [ ] Identify and mark evolved stars in CMD
- [ ] Implement per-star MCMC:
  ```python
  def mcmc_mass_inference(g_obs, g_err, r_obs, r_err, iso, membership_prob)
  ```
- [ ] Add population-level hierarchical model (PyMC3)
- [ ] Implement mass segregation test
- [ ] Create mass histogram with uncertainties
- [ ] Compare isochrones: PARSEC vs. MIST vs. empirical MLR
- [ ] Propagate photometric errors through entire mass pipeline

---

### 5. NEW: data_reduction.ipynb - Status: 0% (CRITICAL - CREATE)

#### Required Core Functions:
```python
# FITS I/O
- read_fits_images(filelist)
- validate_fits_headers(fits_data)
- extract_image_data(fits_hdu)

# Image Preprocessing
- bias_correction(image, bias_frame)
- dark_correction(image, dark_frame, exptime)
- flat_field_correction(image, flat_frame)
- cosmic_ray_rejection(image_list, method='median')
- image_stack_align_combine(image_list, reference_fits, method='median_combine')

# PSF Characterization
- measure_psf(image, star_catalogs, box_size=25)
- fit_psf_model(image, x, y, psf_model='moffat')
- psf_interpolator(psf_measurements, grid)
- apply_variable_psf_correction(image, psf_map)

# Source Detection & Photometry
- detect_sources(image, threshold, fwhm, min_area)
- centroid_sources(image, x, y, method='center_of_mass')
- psf_photometry(image, x, y, psf_model, uncertainties=True)
- aperture_photometry(image, x, y, aperture_radius, background_method='local')
- photometry_quality_assessment(photometry_cat)

# Astrometry
- wcs_solution_gaia(image, gaia_catalog, order=2)
- refine_wcs_solution(image, catalog, method='twangs')
- apply_wcs_to_catalog(photometry_cat, wcs_solution)
- validate_astrometry(catalog, gaia_reference)

# Photometric Calibration
- find_photometric_standards(catalog, sdss_reference)
- photometric_zeropoint(catalog, standard_mags)
- apply_photometric_calibration(catalog, zeropoint, color_term=None)
- magnitude_to_flux_conversion(mag, mag_err, zeropoint)
```

#### Implementation Strategy:
1. Use `astropy.io.fits` for FITS handling
2. Use `photutils.background` for background subtraction
3. Use `astropy.convolution` for PSF modeling
4. Use `photutils.detection` for source finding
5. Use `photutils.psf` for PSF photometry
6. Use `photutils.aperture` for aperture photometry
7. Use `astropy.wcs` for astrometric solutions
8. Use `astroquery.vizier` for SDSS standard calibration

---

### 6. NEW: artificial_star_tests.ipynb - Status: 0% (CRITICAL - CREATE)

#### Required Functions:
```python
def inject_artificial_stars(image, magnitude_bins, positions, n_per_bin):
    """
    Inject artificial stars into image copies.
    Returns: (injected_images, true_catalog)
    """

def run_photometry_on_artificial(injected_images, photometry_func):
    """
    Run photometry pipeline on artificial images.
    Returns: detected_catalog
    """

def measure_recovery_fraction(true_catalog, detected_catalog, mag_bins):
    """
    Compare detected to true -> recovery fraction C(m)
    Returns: C(m), recovery_counts
    """

def spatial_completeness_map(image_shape, photometry_results):
    """
    Create 2D completeness map (how it varies spatially)
    Returns: completeness_map, uncertainties
    """

def color_dependent_completeness(color_bins, catalog):
    """
    Measure if completeness varies with color
    Returns: C(m, color)
    """
```

#### Workflow:
1. Create artificial stars across magnitude range (15-22 mag)
2. Inject into actual science images
3. Run same photometry pipeline
4. Measure recovery rates in magnitude/position/color bins
5. Fit completeness functions (all three models)
6. Compare models via information criteria
7. Test spatial/color dependence

---

### 7. NEW: density_profiles.ipynb - Status: 0% (CRITICAL - CREATE)

#### Required Functions:
```python
# Radial Binning
def radial_bin_stars(x, y, center_x, center_y, bin_edges='auto', method='linear'):
    """Create radial bins and assign stars"""

def background_surface_density(x, y, bin_edges, outer_radius_percentile=90):
    """Estimate background from outer regions"""

def surface_density_profile(membership_cat, completeness_func, bin_edges):
    """Construct Σ(R) with completeness corrections"""

def radial_profile_errors(membership_cat, completeness_cat, bin_edges, method='poisson'):
    """Propagate uncertainties through profile"""

# Model Fitting
def plummer_profile(r, M_total, a):
    """Plummer surface density (Equation 27)"""

def king_profile(r, rho0, rc, rt):
    """King profile surface density"""

def likelihood_profile_fit(bin_centers, surface_density, errors, profile_func, params):
    """Maximum likelihood profile fitting"""

def fit_plummer(radii, surface_density, errors):
    """Fit Plummer model, return (M, a) with uncertainties"""

def fit_king(radii, surface_density, errors):
    """Fit King model, return (rc, rt, c) with uncertainties"""

# Model Comparison
def model_comparison_metrics(data, errors, profile_fits):
    """Compute AIC, BIC, likelihood ratios"""

def bootstrap_profile_uncertainties(membership_cat, completeness_func, 
                                    bin_edges, n_bootstrap=1000):
    """Bootstrap resampling of profile parameters"""

# Visualization
def plot_density_profile(radii, sigma, sigma_err, models):
    """Plot observed profile + model fits + residuals"""
```

#### Analysis Steps:
1. Define radial bins (adjustable spacing)
2. Estimate background from outer field
3. Construct surface density profile with weights
4. Fit Plummer, King, Wilson models
5. Compare models (AIC/BIC/likelihood ratios)
6. Bootstrap for uncertainties
7. Extrapolate to total cluster mass
8. Compare M2 vs. M34 profiles

---

### 8. NEW: bayesian_hierarchical.ipynb - Status: 0% (ADVANCED)

#### Hierarchical Model Structure:
```python
# Cluster parameters (fixed per iteration but updated in hierarchy)
μ_cluster ~ Normal(μ_prior, σ_prior)
Σ_cluster ~ InverseWishart(ν, Ψ)

# Field parameters
μ_field ~ Normal(μ_field_prior, σ_field_prior)
Σ_field ~ InverseWishart(ν, Ψ)

# Population parameters (shared across all stars)
f_member ~ Beta(2, 2)                    # Membership fraction
α_imf ~ Normal(2.3, 0.5)                # IMF slope
f_binary ~ Beta(5, 5)                   # Binary fraction
M_turn_off ~ Normal(0.8, 0.2)          # Main seq turnoff mass

# Per-star membership
z_i ~ Bernoulli(f_member)

# Per-star data likelihood
g_i | z_i ~ { 
    if z_i=1: N(g_expected_iso(M_i), σ_g)    # Cluster
    if z_i=0: Uniform(15, 25)                  # Field
}

pm_i | z_i ~ {
    if z_i=1: N(μ_cluster, Σ_cluster)
    if z_i=0: N(μ_field, Σ_field)
}

r_i | z_i ~ {
    if z_i=1: Plummer(a, M)
    if z_i=0: Uniform(field)
}
```

#### Implementation:
```python
import pymc3 as pm

with pm.Model() as hierarchical_model:
    # Priors on population parameters
    f_member = pm.Beta('f_member', 2, 2)
    
    # Cluster kinematics
    mu_cl = pm.Normal('mu_cluster', mu=0, sigma=10, shape=2)
    sigma_cl = pm.HalfNormal('sigma_cluster', sigma=2, shape=2)
    
    # Membership
    z = pm.Bernoulli('membership', p=f_member, shape=N_stars)
    
    # Observations and likelihood
    # ... (detailed PM, CMD, spatial likelihoods)
    
    trace = pm.sample(2000, tune=1000)
```

---

### 9. NEW: systematic_errors.ipynb - Status: 0% (IMPORTANT)

#### Error Budget Tracking:
```python
class SystematicErrorBudget:
    """Track and propagate systematic uncertainties"""
    
    def __init__(self):
        self.sources = {}
    
    def add_photometric_errors(self, catalog, mag_err, color_err):
        """Photometric calibration uncertainty"""
        
    def add_completeness_errors(self, completeness, completeness_err):
        """Uncertainty in C(m) fit"""
        
    def add_membership_errors(self, membership_cat):
        """Membership probability uncertainties"""
        
    def add_distance_error(self, distance_pc, distance_err):
        """Distance modulus uncertainty"""
        
    def add_extinction_error(self, Av, Av_err):
        """Reddening uncertainty"""
        
    def add_isochrone_error(self, age_err, feh_err):
        """Age/metallicity uncertainty effect"""
        
    def propagate_to_profile(self, profile_func):
        """Monte Carlo error propagation through analysis chain"""
        return profile_profile_err
        
    def generate_budget_report(self):
        """Print summary of error contributions"""
```

#### Propagation Methods:
- Monte Carlo sampling (N=10000)
- Analytical covariance propagation
- Bootstrap resampling
- Jackknife analysis
- Sensitivity analysis

---

## INTEGRATION WORKFLOW

### Dependency Graph:
```
data_reduction.ipynb
    ↓ (produces: photometry_catalog.fits)
    ├→ artificial_star_tests.ipynb → completeness.ipynb
    ├→ membership.ipynb (needs photometry + Gaia)
    └→ mass_estimation.ipynb
    
    ↓
    density_profiles.ipynb
    
    ├→ systematic_errors.ipynb
    ├→ bayesian_hierarchical.ipynb
    └→ (output: final results)
```

### Master Orchestration:
```python
# analysis_pipeline.ipynb
from importlib import reload

# Configuration
config = {
    'cluster_name': 'M2',
    'cluster_coords': (323.3626, -0.8233),
    'age_gyr': 13.0,
    'metallicity_feh': -1.6,
    'distance_pc': 11500,
    'av': 0.186,
    'psf_photometry': True,
}

# Execute pipeline
photometry = data_reduction.process_images(config)
completeness = artificial_star_tests.fit_completeness(photometry, config)
membership = membership.determine_membership(photometry, completeness, config)
masses = mass_estimation.estimate_masses(photometry, membership, config)
profiles = density_profiles.construct_profiles(membership, completeness, config)
errors = systematic_errors.propagate_errors(photometry, membership, profiles, config)

# Results
results = {
    'photometry': photometry,
    'completeness': completeness,
    'membership': membership,
    'masses': masses,
    'profiles': profiles,
    'uncertainties': errors
}
```

---

## TESTING & VALIDATION STRATEGY

### Unit Tests:
- [ ] Test each function with synthetic data
- [ ] Test edge cases (empty bins, single stars, etc.)
- [ ] Test error propagation with known uncertainties

### Integration Tests:
- [ ] Run full pipeline on mock cluster data
- [ ] Compare to published results for known clusters
- [ ] Cross-validate with independent methods

### Validation Against Known Data:
- [ ] Compare membership to published cluster catalogs
- [ ] Compare profiles to GAIA DR3 density profiles
- [ ] Compare masses to spectroscopic masses (if available)

---

## DOCUMENTATION REQUIREMENTS

### Per Notebook:
- [ ] Function docstrings (numpy format)
- [ ] Algorithm descriptions with references
- [ ] Example usage cells
- [ ] Validation plots
- [ ] Assumptions and limitations section

### Project-Level:
- [ ] README with overview
- [ ] METHODOLOGY.md with detailed methodology
- [ ] DATA.md describing input/output formats
- [ ] REFERENCES.md with all citations
- [ ] TROUBLESHOOTING.md with common issues

---

## ESTIMATED TIME ALLOCATION

| Task | Duration | Priority |
|------|----------|----------|
| Data Reduction | 2 weeks | CRITICAL |
| Artificial Star Tests | 1 week | CRITICAL |
| Real Data Validation | 1 week | CRITICAL |
| Density Profiles | 1 week | CRITICAL |
| Mass Estimation (real isochrones) | 1 week | IMPORTANT |
| Membership Refinement | 1 week | IMPORTANT |
| Systematic Errors | 1 week | IMPORTANT |
| Bayesian Hierarchical | 1 week | ADVANCED |
| Integration & Testing | 1 week | CRITICAL |
| Documentation | 1 week | IMPORTANT |
| Paper Writing | 2 weeks | CRITICAL |
| **TOTAL** | **14 weeks** | - |

---

## SUCCESS CRITERIA

### Minimum Viable Product (MVP):
✓ Complete photometric catalog (M2 + M34)
✓ Verified completeness functions
✓ Membership probabilities on real data
✓ Density profiles with model fits
✓ Publication-ready figures

### Full Implementation:
✓ All gaps closed
✓ Full systematic error budget
✓ Advanced statistical methods
✓ Comprehensive documentation
✓ Submission-ready manuscript

### Excellence Tier:
✓ Hierarchical Bayesian modeling
✓ Mass segregation analysis
✓ Machine learning validation
✓ Dynamical modeling
✓ Software package release

