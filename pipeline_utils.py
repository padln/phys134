"""
Data Reduction Pipeline Utilities for M34 Photometry

This module provides clean, reusable functions for:
- Background estimation
- Source detection
- Aperture photometry
- Photometric calibration
- Catalog generation

Usage in notebook:
    from pipeline_utils import *
    data_sub, bkg = estimate_background(data, box_size=64)
    sources = detect_sources(data_sub, bkg, fwhm=3.5, threshold=3.0)
"""

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import SigmaClip
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from photutils.background import Background2D, MedianBackground
from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry


def load_fits_image(filepath):
    """
    Load FITS file and extract data, header, and WCS.

    Parameters:
    -----------
    filepath : str
        Path to FITS file

    Returns:
    --------
    data : ndarray
        Image data (2D array)
    header : Header
        FITS header
    wcs : WCS
        World Coordinate System object
    """
    with fits.open(filepath) as hdul:
        data = hdul[0].data.astype(np.float64)
        header = hdul[0].header
        wcs = WCS(header)

    return data, header, wcs


def estimate_background(data, box_size=64, filter_size=3, sigma=3.0, maxiters=5):
    """
    Estimate and subtract 2D background with sigma clipping.

    Parameters:
    -----------
    data : ndarray
        Image data
    box_size : int
        Size of background mesh boxes (pixels)
    filter_size : int
        Median filter size for smoothing
    sigma : float
        Sigma clipping threshold
    maxiters : int
        Maximum sigma clipping iterations

    Returns:
    --------
    data_sub : ndarray
        Background-subtracted image
    bkg : Background2D
        Background object with .background and .background_rms
    """
    sigma_clip = SigmaClip(sigma=sigma, maxiters=maxiters)
    bkg_estimator = MedianBackground()

    bkg = Background2D(data, box_size=box_size, filter_size=filter_size,
                       sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)

    data_sub = data - bkg.background

    return data_sub, bkg


def detect_sources(data_sub, bkg, fwhm, threshold=3.0, sharplo=0.2, sharphi=1.0,
                   roundlo=-1.0, roundhi=1.0):
    """
    Detect point sources using DAOStarFinder.

    Parameters:
    -----------
    data_sub : ndarray
        Background-subtracted image
    bkg : Background2D
        Background object (for RMS)
    fwhm : float
        FWHM of PSF in pixels
    threshold : float
        Detection threshold in units of background RMS
    sharplo, sharphi : float
        Sharpness limits (reject cosmic rays)
    roundlo, roundhi : float
        Roundness limits (reject elongated objects)

    Returns:
    --------
    sources : Table
        Detected sources with positions, fluxes, etc.
    """
    daofind = DAOStarFinder(fwhm=fwhm,
                            threshold=threshold * bkg.background_rms_median,
                            sharplo=sharplo, sharphi=sharphi,
                            roundlo=roundlo, roundhi=roundhi)

    sources = daofind(data_sub)

    if sources is not None:
        # Add instrumental magnitudes
        sources['mag_inst'] = -2.5 * np.log10(sources['flux'])
        # Sort by brightness
        sources.sort('flux', reverse=True)

    return sources


def aperture_photometry_pipeline(data_sub, sources, aperture_radius,
                                  annulus_inner, annulus_outer,
                                  gain, rdnoise, exptime):
    """
    Perform aperture photometry with proper uncertainty propagation.

    Parameters:
    -----------
    data_sub : ndarray
        Background-subtracted image
    sources : Table
        Source positions (must have 'xcentroid', 'ycentroid')
    aperture_radius : float
        Aperture radius in pixels
    annulus_inner, annulus_outer : float
        Sky annulus radii in pixels
    gain : float
        Detector gain (e-/ADU)
    rdnoise : float
        Read noise (electrons)
    exptime : float
        Exposure time (seconds)

    Returns:
    --------
    phot_table : Table
        Photometry results with fluxes and uncertainties
    """
    positions = np.transpose([sources['xcentroid'], sources['ycentroid']])

    # Define apertures
    aperture = CircularAperture(positions, r=aperture_radius)
    annulus = CircularAnnulus(positions, r_in=annulus_inner, r_out=annulus_outer)

    # Measure flux in apertures and annuli
    phot_table = aperture_photometry(data_sub, [aperture, annulus])

    # Calculate mean sky background per pixel from annulus
    annulus_area = annulus.area
    sky_per_pixel = phot_table['aperture_sum_1'] / annulus_area

    # Sky-subtracted flux
    aperture_area = aperture.area
    total_sky = sky_per_pixel * aperture_area
    flux_aper = phot_table['aperture_sum_0'] - total_sky

    # Uncertainty calculation
    # Convert to electrons for noise calculation
    flux_electrons = flux_aper * gain

    # Poisson noise from source
    var_source = flux_electrons

    # Sky background noise (Poisson)
    sky_electrons_per_pix = sky_per_pixel * gain
    var_sky = aperture_area * sky_electrons_per_pix

    # Read noise
    var_read = aperture_area * rdnoise**2

    # Total variance (in electrons)
    var_total_electrons = var_source + var_sky + var_read

    # Convert back to ADU
    var_total_adu = var_total_electrons / gain**2
    flux_err = np.sqrt(var_total_adu)

    # Add to table
    phot_table['flux_aper'] = flux_aper
    phot_table['flux_err'] = flux_err

    # Convert to magnitudes
    # Handle negative/zero fluxes
    good_flux = flux_aper > 0
    mag = np.full(len(flux_aper), np.nan)
    mag_err = np.full(len(flux_aper), np.nan)

    mag[good_flux] = -2.5 * np.log10(flux_aper[good_flux])
    mag_err[good_flux] = (2.5 / np.log(10)) * (flux_err[good_flux] / flux_aper[good_flux])

    phot_table['mag_inst'] = mag
    phot_table['mag_err'] = mag_err

    return phot_table


def cross_match_catalogs(sources_1, sources_2, wcs_1, wcs_2, max_sep_arcsec=2.0):
    """
    Cross-match two source catalogs using sky coordinates.

    Parameters:
    -----------
    sources_1, sources_2 : Table
        Source catalogs with 'xcentroid', 'ycentroid'
    wcs_1, wcs_2 : WCS
        World coordinate systems
    max_sep_arcsec : float
        Maximum separation for match (arcseconds)

    Returns:
    --------
    matches : Table
        Matched sources from both catalogs
    idx_1, idx_2 : ndarray
        Indices of matches in original catalogs
    sep : Quantity
        Separations of matches
    """
    from astropy.coordinates import match_coordinates_sky

    # Convert pixel coordinates to sky coordinates
    coords_1 = wcs_1.pixel_to_world(sources_1['xcentroid'], sources_1['ycentroid'])
    coords_2 = wcs_2.pixel_to_world(sources_2['xcentroid'], sources_2['ycentroid'])

    # Match catalogs
    idx, sep, _ = match_coordinates_sky(coords_1, coords_2)

    # Filter by separation
    good_match = sep < max_sep_arcsec * u.arcsec

    idx_1 = np.where(good_match)[0]
    idx_2 = idx[good_match]
    sep_matched = sep[good_match]

    return idx_1, idx_2, sep_matched


def calculate_photometric_zeropoint(sources, ref_mags, ref_mag_errs, sigma_clip=3.0):
    """
    Calculate photometric zeropoint using reference stars.

    ZP is defined such that: mag_calibrated = mag_inst + ZP

    Parameters:
    -----------
    sources : Table
        Source catalog with 'mag_inst', 'mag_err'
    ref_mags : array
        Reference magnitudes (calibrated)
    ref_mag_errs : array
        Reference magnitude uncertainties
    sigma_clip : float
        Sigma clipping threshold

    Returns:
    --------
    zeropoint : float
        Photometric zeropoint
    zp_err : float
        Uncertainty in zeropoint
    n_used : int
        Number of stars used
    """
    # Calculate zeropoint for each star
    zp_individual = ref_mags - sources['mag_inst']

    # Sigma clip outliers
    from astropy.stats import sigma_clipped_stats
    zp_mean, zp_median, zp_std = sigma_clipped_stats(zp_individual, sigma=sigma_clip)

    # Use median as robust estimator
    zeropoint = zp_median

    # Uncertainty: std / sqrt(N) for clipped sample
    mask = np.abs(zp_individual - zp_median) < sigma_clip * zp_std
    n_used = np.sum(mask)
    zp_err = zp_std / np.sqrt(n_used) if n_used > 0 else np.inf

    return zeropoint, zp_err, n_used


def save_catalog(catalog, filename, format='csv', overwrite=True):
    """
    Save source catalog to file.

    Parameters:
    -----------
    catalog : Table
        Source catalog
    filename : str
        Output filename
    format : str
        File format ('csv', 'fits', 'ascii')
    overwrite : bool
        Overwrite existing file
    """
    catalog.write(filename, format=format, overwrite=overwrite)
    print(f"✓ Catalog saved to {filename}")


# Print module info when imported
print("✓ Pipeline utilities loaded")
print("  Available functions:")
print("    - load_fits_image(filepath)")
print("    - estimate_background(data, box_size=64)")
print("    - detect_sources(data_sub, bkg, fwhm, threshold=3.0)")
print("    - aperture_photometry_pipeline(...)")
print("    - cross_match_catalogs(sources_1, sources_2, wcs_1, wcs_2)")
print("    - calculate_photometric_zeropoint(...)")
print("    - save_catalog(catalog, filename)")
