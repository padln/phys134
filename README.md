# Stellar Cluster Density Profile Analysis
## PHYS 134L - Advanced Observational Astrophysics

[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![AstroPy](https://img.shields.io/badge/powered%20by-AstroPy-orange.svg)](http://www.astropy.org/)

A comprehensive photometric analysis pipeline for measuring stellar density profiles in globular and open clusters, with applications to M2 (NGC 7089) and M34 (NGC 1039).

## Overview

This project implements a statistically rigorous pipeline for:
- **CCD photometry** with proper background estimation and uncertainty propagation
- **Completeness corrections** via artificial star tests
- **Bayesian membership determination** combining CMD, proper motion, and spatial information
- **Density profile fitting** using Plummer and King models
- **Observation planning** via Fisher information analysis

## Quick Start

### Prerequisites

```bash
# Python 3.8+ with standard scientific packages
pip install numpy scipy matplotlib astropy photutils astroquery emcee
```

### Running the Pipeline

1. **Data Reduction** - Process raw FITS images to photometric catalog:
   ```bash
   jupyter notebook data_reduction_simple.ipynb
   ```
   - Inputs: g and r-band FITS files in `Data/`
   - Output: `m34_photometry_instrumental.csv`

2. **Membership Determination** - Filter field stars using multi-dimensional criteria:
   ```bash
   jupyter notebook membership.ipynb
   ```
   - Inputs: Photometric catalog, Gaia proper motions
   - Output: Cluster member catalog with probabilities

3. **Density Profile** - Fit Plummer/King models to radial distribution:
   ```bash
   jupyter notebook density_profile.ipynb
   ```
   - Input: Member catalog
   - Output: Best-fit profile parameters, plots

## Repository Structure

```
PHYS 134L/
├── Data/                          # FITS images (not tracked in git)
├── AASTeX_Template/               # LaTeX article + figures
│   └── article.tex                # Research paper
├── docs/                          # Theory and methodology guides
│   ├── BACKGROUND_ESTIMATION_THEORY.md
│   └── PIPELINE_WORKFLOW.md
├── pipeline_utils.py              # Core photometry functions
├── data_reduction_simple.ipynb    # Step 1: FITS → catalog
├── membership.ipynb               # Step 2: Catalog → members
├── completeness.ipynb             # Completeness function analysis
├── signal_noise.ipynb             # S/N calculations
└── mass_estimation.ipynb          # Photometric mass estimates
```

## Key Features

### 🔭 Observation Planning
- Signal-to-noise calculations for exposure time optimization
- Multi-band imaging strategy (SDSS g', r', i')
- Las Cumbres Observatory 0.4m telescope specifications

### 📊 Photometry Pipeline
- 2D background estimation with sigma clipping
- DAOStarFinder source detection
- Aperture photometry with full uncertainty propagation
- Cross-matching between filters

### ⭐ Completeness Modeling
- Artificial star test framework
- Error function, tanh, and Fermi-Dirac model comparison
- Richardson-Lucy deconvolution for luminosity function recovery
- Bayesian model selection (AIC/BIC)

### 🎯 Membership Determination
- **Color-Magnitude Diagram** filtering with isochrone matching
- **Proper Motion** filtering using Gaia DR3 astrometry
- **Spatial distribution** modeling
- Bayesian probability combination

### 📈 Profile Fitting
- Plummer profile: $\rho(r) = \rho_0 (1 + r^2/a^2)^{-5/2}$
- King profile for comparison
- Maximum likelihood and MCMC parameter estimation

## Theoretical Background

See [`docs/`](docs/) for detailed theoretical explanations:
- **Background Estimation**: Algorithm details, validation, impact on photometry
- **Pipeline Workflow**: How all the pieces fit together

## Example Results

### M34 Color-Magnitude Diagram
Proper motion-selected cluster members showing clear main sequence:

![CMD placeholder]

### Background Estimation
2D mesh-based background model capturing spatial gradients:

![Background panel placeholder - see AASTeX_Template/background_panel.png]

### Completeness Function
Three functional forms showing excellent agreement:

![Completeness comparison placeholder]

## Data Sources

- **Imaging**: Las Cumbres Observatory 0.4m network (SDSS g', r' filters)
- **Astrometry**: Gaia DR3 (proper motions, parallaxes)
- **Photometric calibration**: Gaia DR3 / Pan-STARRS DR2

## Scientific Context

This project compares two stellar systems with vastly different properties:

| Property | M2 (Globular) | M34 (Open) |
|----------|---------------|------------|
| Age | ~13 Gyr | ~200 Myr |
| Mass | ~10⁵ M☉ | ~10³ M☉ |
| Distance | 11.5 kpc | 470 pc |
| [Fe/H] | -1.6 dex | 0.0 dex |
| Concentration | High (c=1.6) | Low (loose) |

The contrasting density profiles reveal how stellar systems evolve under different dynamical conditions.

## Publication

Full methodology and results are described in:
- **Paper**: `AASTeX_Template/article.tex` (compile with `pdflatex` + `bibtex`)

## Citation

If you use this code for your research, please cite:
```bibtex
@misc{phys134l_cluster_profiles,
  author = {[Nathan Madsen]},
  title = {Stellar Cluster Density Profile Analysis},
  year = {2025},
  publisher = {GitHub},
  url = {https://github.com/padln/phys134}
}
```

## Dependencies

Core packages:
- `numpy` - Numerical computing
- `scipy` - Optimization and statistics
- `matplotlib` - Visualization
- `astropy` - FITS I/O, WCS, units, coordinates
- `photutils` - Background estimation, aperture photometry
- `astroquery` - Gaia data access
- `emcee` - MCMC sampling

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file for details.

## Acknowledgments

- **Las Cumbres Observatory** for telescope time
- **Gaia mission** for astrometric data
- **UCSB Physics Department** for course support
- **Photutils development team** for excellent documentation

## Contact

For questions about this code or methodology:
- Open an issue on GitHub
- Email: [madsen@ucsb.edu]

---

**Course**: PHYS 134L - Advanced Observational Astrophysics
**Institution**: UC Santa Barbara
**Term**: Fall 2025
