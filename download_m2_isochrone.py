#!/usr/bin/env python3
"""
Download PARSEC isochrone for M2 globular cluster

M2 (NGC 7089) parameters:
- Age: 13 Gyr (logAge = 10.11)
- [Fe/H] = -1.6 (Z ≈ 0.0001)
- Distance: 11.5 kpc (distance modulus μ = 15.30 mag)
- E(B-V) = 0.06 mag
"""

print("=" * 70)
print("M2 ISOCHRONE DOWNLOAD INSTRUCTIONS")
print("=" * 70)
print("\nGo to: http://stev.oapd.inaf.it/cgi-bin/cmd")
print("\nSettings:")
print("  1. Evolutionary tracks: PARSEC version 1.2S + COLIBRI")
print("  2. Photometric system: SDSS ugriz (Vega or AB)")
print("  3. Circumstellar dust: OFF")
print("  4. Interstellar extinction: Av = 0 (we'll apply manually)")
print("\n  5. Single isochrone with:")
print("     - log(age/yr) = 10.11  (13 Gyr)")
print("     - [M/H] = -1.6  (or Z = 0.0001 if using Z)")
print("     - Initial mass range: 0.1 - 10 Msun")
print("\n  6. Output: Long format, no dust")
print("\n  7. Download and save as: isochrone_m2_13gyr_mh-1.6.dat")
print("\n  8. Place in current directory")
print("=" * 70)

print("\nAlternatively, try this direct URL (may need adjustment):")
print("http://stev.oapd.inaf.it/cgi-bin/cmd_3.7")

print("\nAfter downloading, we'll apply:")
print(f"  - Distance modulus: μ = 15.30 mag (11.5 kpc)")
print(f"  - Reddening: E(B-V) = 0.06 mag")
print(f"  - A_g = 0.20 mag, A_r = 0.14 mag")
print("=" * 70)
