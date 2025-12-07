#!/usr/bin/env python3
"""
Run Artificial Star Test on M2 and M34 data
"""
import numpy as np
import warnings
warnings.filterwarnings('ignore')

from pipeline_utils import (
    load_fits_image,
    run_artificial_star_test,
    print_ast_results
)
from scipy.optimize import minimize
from scipy.special import erf
import os

def completeness_model(m, m50, sigma_comp):
    """Error function completeness model"""
    return 0.5 * (1 + erf((m50 - m) / (np.sqrt(2) * sigma_comp)))

def negative_log_likelihood(theta, m_bins, N_add, N_rec):
    """Negative log-likelihood for binomial model"""
    m50, sigma_comp = theta
    C = completeness_model(m_bins, m50, sigma_comp)
    C = np.clip(C, 1e-10, 1 - 1e-10)
    log_L = np.sum(N_rec * np.log(C) + (N_add - N_rec) * np.log(1 - C))
    return -log_L

def run_ast_for_cluster(fits_file, cluster_name, output_prefix):
    """Run AST for a single cluster image"""
    print(f"\n{'='*70}")
    print(f"Running AST for {cluster_name}")
    print(f"File: {fits_file}")
    print(f"{'='*70}\n")
    
    # Load image
    data, header, wcs = load_fits_image(fits_file)
    print(f"✓ Loaded image: {data.shape}")
    print(f"  Exposure time: {header.get('EXPTIME', 'N/A')} s")
    print(f"  Filter: {header.get('FILTER', 'N/A')}")
    
    # Define magnitude bins - use realistic range for instrumental mags
    # Based on analysis: real stars are around -15 to -8 mag instrumental
    mag_bins = np.arange(-16.0, -7.0, 1.0)
    n_stars_per_bin = 100
    fwhm = header.get('L1FWHM', 3.5)
    
    print(f"\nAST Configuration:")
    print(f"  Magnitude range: {mag_bins.min():.1f} to {mag_bins.max():.1f}")
    print(f"  Number of bins: {len(mag_bins)}")
    print(f"  Stars per bin: {n_stars_per_bin}")
    print(f"  Total stars: {len(mag_bins) * n_stars_per_bin}")
    print(f"  FWHM: {fwhm:.2f} pixels")
    print(f"\nRunning AST... (this may take a few minutes)")
    
    # Run AST
    results = run_artificial_star_test(
        data=data,
        header=header,
        wcs=wcs,
        magnitude_bins=mag_bins,
        n_stars_per_bin=n_stars_per_bin,
        fwhm=fwhm,
        zeropoint=0.0,
        match_radius=2.0
    )
    
    # Print results
    print_ast_results(results)
    
    # Fit completeness model
    m_bins = results['magnitude_bins']
    N_add = results['n_added']
    N_rec = results['n_recovered']
    C_obs = results['completeness']
    C_err = results['completeness_err']
    
    # Initial guess: m50 in middle of range, sigma ~ 1 mag
    initial_guess = [np.mean(mag_bins), 1.0]
    
    result_fit = minimize(
        negative_log_likelihood,
        initial_guess,
        args=(m_bins, N_add, N_rec),
        method='Nelder-Mead'
    )
    
    m50_fit, sigma_fit = result_fit.x
    
    print(f"\n{'='*70}")
    print("Fitted Completeness Model:")
    print(f"{'='*70}")
    print(f"  m50 (50% completeness): {m50_fit:.2f} mag")
    print(f"  sigma_comp (width): {sigma_fit:.2f} mag")
    print(f"  Log-likelihood: {-result_fit.fun:.2f}")
    
    # Calculate some useful metrics
    # The value 0.906 comes from solving erf^-1(0.8) / sqrt(2)
    # for 90% completeness: C = 0.5 * (1 + erf((m50 - m90) / (sqrt(2) * sigma)))
    # Solving for m90 when C=0.9 gives m90 = m50 + sigma * sqrt(2) * 0.906
    ERF_INV_08 = 0.906  # erf^-1(0.8) / sqrt(2), for 90%/10% completeness levels
    m90 = m50_fit + sigma_fit * np.sqrt(2) * ERF_INV_08
    m10 = m50_fit - sigma_fit * np.sqrt(2) * ERF_INV_08
    
    print(f"\nCompleteness Levels:")
    print(f"  90% complete at: {m90:.2f} mag")
    print(f"  50% complete at: {m50_fit:.2f} mag")
    print(f"  10% complete at: {m10:.2f} mag")
    print(f"  Transition width: {m10 - m90:.2f} mag")
    print(f"{'='*70}\n")
    
    # Save results
    results_file = f"{output_prefix}_ast_results.txt"
    params_file = f"{output_prefix}_model_params.txt"
    
    with open(results_file, 'w') as f:
        f.write(f"# Artificial Star Test Results - {cluster_name}\n")
        f.write(f"# Image: {fits_file}\n")
        f.write("# Fitted Completeness Model: Error Function\n")
        f.write(f"# m50 = {m50_fit:.4f} mag\n")
        f.write(f"# sigma_comp = {sigma_fit:.4f} mag\n")
        f.write("#\n")
        f.write("# Magnitude  N_added  N_recovered  Completeness  Error\n")
        
        for i in range(len(m_bins)):
            f.write(f"{m_bins[i]:10.3f}  {N_add[i]:7d}  {N_rec[i]:11d}  "
                   f"{C_obs[i]:12.4f}  {C_err[i]:7.4f}\n")
    
    print(f"✓ Results saved to {results_file}")
    
    with open(params_file, 'w') as f:
        f.write(f"# Completeness model parameters for {cluster_name}\n")
        f.write(f"m50 {m50_fit:.4f}\n")
        f.write(f"sigma_comp {sigma_fit:.4f}\n")
    
    print(f"✓ Model parameters saved to {params_file}")
    
    return results, m50_fit, sigma_fit

# Main execution
if __name__ == "__main__":
    data_dir = "./Data"
    
    # Find FITS files
    fits_files = []
    if os.path.exists(data_dir):
        for root, dirs, files in os.walk(data_dir):
            for f in files:
                if f.endswith(".fits") and not f.endswith(".fits.fz"):
                    fits_files.append(os.path.join(root, f))
    
    if len(fits_files) == 0:
        print("No FITS files found in Data/ directory")
        exit(1)
    
    print(f"Found {len(fits_files)} FITS file(s):")
    for i, f in enumerate(fits_files):
        print(f"  {i+1}. {f}")
    
    # Run AST on first two files
    # Try to identify cluster from FITS header or filename
    for i, fits_file in enumerate(fits_files[:2]):
        # Try to extract cluster name from header
        try:
            from astropy.io import fits as pyfits
            with pyfits.open(fits_file) as hdul:
                header = hdul[0].header
                obj_name = header.get('OBJECT', '').upper()
                filter_name = header.get('FILTER', 'unknown')
                
                # Identify cluster from OBJECT header
                if 'M2' in obj_name or 'NGC 7089' in obj_name or 'NGC7089' in obj_name:
                    cluster_name = 'M2'
                elif 'M34' in obj_name or 'NGC 1039' in obj_name or 'NGC1039' in obj_name:
                    cluster_name = 'M34'
                else:
                    # Fall back to generic naming
                    cluster_name = f'Cluster_{i+1}'
                
                # Add filter info to output prefix
                output_prefix = f"{cluster_name.lower()}_{filter_name}_completeness"
        except Exception as e:
            print(f"Warning: Could not read header from {fits_file}: {e}")
            cluster_name = f'Cluster_{i+1}'
            output_prefix = f"{cluster_name.lower()}_completeness"
        
        run_ast_for_cluster(fits_file, cluster_name, output_prefix)
    
    print("\n" + "="*70)
    print("AST runs complete!")
    print("="*70)
