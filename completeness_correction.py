"""
Completeness Correction Module for Radial Profiles

This module provides functions to:
1. Load completeness model parameters from AST results
2. Apply completeness corrections to observed star counts
3. Apply Richardson-Lucy deconvolution for improved profile recovery
"""

import numpy as np
from scipy.special import erf
from scipy.ndimage import convolve1d


def load_completeness_model(params_file='completeness_model_params.txt'):
    """
    Load completeness model parameters from file.
    
    Parameters:
    -----------
    params_file : str
        Path to file containing m50 and sigma_comp parameters
        
    Returns:
    --------
    m50 : float
        Magnitude at 50% completeness
    sigma_comp : float
        Completeness transition width
    """
    params = {}
    with open(params_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split()
            if len(parts) == 2:
                params[parts[0]] = float(parts[1])
    
    return params.get('m50'), params.get('sigma_comp')


def completeness_function(m, m50, sigma_comp):
    """
    Error function completeness model.
    
    C(m) = 0.5 * (1 + erf((m50 - m) / (sqrt(2) * sigma_comp)))
    
    Parameters:
    -----------
    m : array-like
        Magnitudes at which to evaluate completeness
    m50 : float
        Magnitude at 50% completeness
    sigma_comp : float
        Completeness transition width
        
    Returns:
    --------
    completeness : array
        Completeness values between 0 and 1
    """
    return 0.5 * (1 + erf((m50 - m) / (np.sqrt(2) * sigma_comp)))


def apply_completeness_correction(N_obs, magnitudes, m50, sigma_comp, min_completeness=0.1):
    """
    Apply completeness correction to observed star counts.
    
    N_true(m) = N_obs(m) / C(m)
    
    Parameters:
    -----------
    N_obs : array
        Observed number of stars in each bin
    magnitudes : array
        Mean magnitude in each bin
    m50 : float
        Magnitude at 50% completeness
    sigma_comp : float
        Completeness transition width
    min_completeness : float
        Minimum completeness to apply correction (avoids division by very small numbers)
        
    Returns:
    --------
    N_corrected : array
        Completeness-corrected star counts
    C : array
        Completeness values used for correction
    """
    # Calculate completeness at each magnitude
    C = completeness_function(magnitudes, m50, sigma_comp)
    
    # Avoid dividing by very small completeness values
    C_safe = np.maximum(C, min_completeness)
    
    # Apply correction
    N_corrected = N_obs / C_safe
    
    return N_corrected, C


def richardson_lucy_deconvolution(observed_profile, psf, n_iterations=10):
    """
    Apply Richardson-Lucy deconvolution to recover true radial profile.
    
    The observed profile is a convolution of the true profile with the PSF
    (point spread function), which includes effects like:
    - Incompleteness as a function of radius (crowding)
    - Spatial variations in detection efficiency
    
    Parameters:
    -----------
    observed_profile : array
        Observed surface density profile
    psf : array
        Point spread function (completeness kernel)
        Should have same length as observed_profile
    n_iterations : int
        Number of Richardson-Lucy iterations
        
    Returns:
    --------
    deconvolved_profile : array
        Recovered true profile
    """
    # Initialize with observed profile
    estimate = observed_profile.copy()
    
    # Normalize PSF
    psf_norm = psf / np.sum(psf)
    
    for i in range(n_iterations):
        # Convolve current estimate with PSF
        convolved = convolve1d(estimate, psf_norm, mode='constant', cval=0.0)
        
        # Avoid division by zero
        convolved_safe = np.maximum(convolved, 1e-10)
        
        # Compute ratio of observed to convolved
        ratio = observed_profile / convolved_safe
        
        # Convolve ratio with flipped PSF
        correction = convolve1d(ratio, psf_norm[::-1], mode='constant', cval=0.0)
        
        # Update estimate
        estimate = estimate * correction
        
        # Ensure non-negativity
        estimate = np.maximum(estimate, 0)
    
    return estimate


def create_completeness_kernel(r_bins, m50, sigma_comp, mag_at_radius_func,
                                fwhm_spatial=1.0):
    """
    Create a spatial completeness kernel for Richardson-Lucy deconvolution.
    
    This accounts for how completeness varies with radius due to crowding
    and other spatial effects.
    
    Parameters:
    -----------
    r_bins : array
        Radial bin edges (in same units as radial profile)
    m50 : float
        Magnitude at 50% completeness
    sigma_comp : float
        Completeness transition width
    mag_at_radius_func : callable
        Function that returns typical magnitude at each radius
        mag = mag_at_radius_func(r)
    fwhm_spatial : float
        FWHM of spatial smoothing kernel (in units of r_bins)
        
    Returns:
    --------
    kernel : array
        Completeness kernel for deconvolution
    """
    r_centers = 0.5 * (r_bins[1:] + r_bins[:-1])
    
    # Get typical magnitude at each radius
    mags = mag_at_radius_func(r_centers)
    
    # Calculate completeness at each radius
    C = completeness_function(mags, m50, sigma_comp)
    
    # Apply spatial smoothing (Gaussian kernel)
    sigma_spatial = fwhm_spatial / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    x = np.arange(len(C))
    kernel = np.exp(-(x - len(C)/2)**2 / (2 * sigma_spatial**2))
    
    # Convolve completeness with spatial kernel
    psf = convolve1d(C, kernel, mode='constant', cval=C[-1])
    
    return psf


def apply_full_correction(r_bins, N_obs, magnitudes, m50, sigma_comp,
                          use_richardson_lucy=True, n_iterations=10,
                          mag_at_radius_func=None):
    """
    Apply full completeness correction including Richardson-Lucy deconvolution.
    
    Parameters:
    -----------
    r_bins : array
        Radial bin edges
    N_obs : array
        Observed star counts in radial bins
    magnitudes : array
        Mean magnitude in each radial bin
    m50 : float
        Magnitude at 50% completeness
    sigma_comp : float
        Completeness transition width
    use_richardson_lucy : bool
        Whether to apply Richardson-Lucy deconvolution
    n_iterations : int
        Number of RL iterations
    mag_at_radius_func : callable, optional
        Function mag = f(r) for creating spatial kernel
        If None, uses provided magnitudes directly
        
    Returns:
    --------
    N_corrected : array
        Completeness-corrected star counts
    N_deconvolved : array or None
        Richardson-Lucy deconvolved counts (if use_richardson_lucy=True)
    """
    # First apply simple completeness correction
    N_corrected, C = apply_completeness_correction(N_obs, magnitudes, m50, sigma_comp)
    
    # Optionally apply Richardson-Lucy deconvolution
    N_deconvolved = None
    if use_richardson_lucy and mag_at_radius_func is not None:
        # Create spatial completeness kernel
        psf = create_completeness_kernel(r_bins, m50, sigma_comp, 
                                         mag_at_radius_func)
        
        # Apply deconvolution to observed counts (before simple correction)
        N_deconvolved = richardson_lucy_deconvolution(N_obs, psf, n_iterations)
    
    return N_corrected, N_deconvolved


# Example usage
if __name__ == "__main__":
    print("Completeness Correction Module")
    print("=" * 70)
    
    # Load completeness model
    try:
        m50, sigma_comp = load_completeness_model()
        print(f"\n✓ Loaded completeness model:")
        print(f"  m50 = {m50:.2f} mag")
        print(f"  sigma_comp = {sigma_comp:.2f} mag")
    except Exception as e:
        print(f"\n✗ Could not load completeness model: {e}")
        print("  Using default values for demonstration")
        m50, sigma_comp = -12.0, 1.0
    
    # Example: correct some radial profile data
    print(f"\n--- Example Correction ---")
    
    # Simulate observed radial profile
    r_bins = np.logspace(0, 2, 11)  # 0.1 to 100 arcmin
    r_centers = 0.5 * (r_bins[1:] + r_bins[:-1])
    
    # Observed counts (decreasing with radius)
    N_obs = 100 * np.exp(-r_centers / 20)
    
    # Typical magnitudes (fainter with radius due to crowding/depth)
    mags = -15 + 0.5 * np.log10(r_centers)
    
    # Apply correction
    N_corr, C = apply_completeness_correction(N_obs, mags, m50, sigma_comp)
    
    print(f"\n{'Radius':>10} | {'N_obs':>10} | {'Mag':>8} | {'C':>8} | {'N_corr':>10}")
    print("-" * 70)
    for i in range(len(r_centers)):
        print(f"{r_centers[i]:10.2f} | {N_obs[i]:10.1f} | {mags[i]:8.2f} | "
              f"{C[i]:8.3f} | {N_corr[i]:10.1f}")
    
    print("\n✓ Completeness correction module ready for use!")
    print("\nTo use in membership analysis:")
    print("  from completeness_correction import apply_completeness_correction")
    print("  N_corrected, C = apply_completeness_correction(N_obs, mags, m50, sigma_comp)")
