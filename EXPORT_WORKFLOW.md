# Export Workflow: membership.ipynb → mass_estimation.ipynb

## Quick Start: Add This Cell at End of membership.ipynb

```python
# ============================================================================
# EXPORT CORRECTED RESULTS FOR MASS_ESTIMATION.IPYNB
# ============================================================================

import os

print("\n" + "="*70)
print("EXPORTING CORRECTED CLUSTER PARAMETERS")
print("="*70)

# Create export dictionary with all relevant quantities
export_data = {
    'cluster_name': cluster_name,
    'age_gyr': age_gyr,
    'distance_pc': distance_pc,
    'parallax_mas': 1000.0 / distance_pc,
    
    # Corrected member counts
    'total_members_uncorrected': int(np.sum(P_combined_safe)),
    'total_members_corrected': int(np.sum(weights_corrected)),
    'correction_factor': correction_factor,
    
    # Radial profile
    'r_centers': r_centers,
    'sigma_members_corrected': sigma_members_corrected,
    'sigma_members_uncorrected': sigma_members_uncorrected,
    'sigma_err_corrected': sigma_err_corrected,
    
    # Luminosity function
    'mag_centers': mag_centers,
    'N_lf_corrected': N_lf_corrected,
    'N_lf_uncorrected': N_lf_uncorrected,
    'N_lf_err_corrected': N_lf_err_corrected,
    
    # IMF
    'mass_centers_imf': mass_centers_imf,
    'N_imf_corrected': N_imf_corrected,
    'N_imf_uncorrected': N_imf_uncorrected,
    'N_imf_err_corrected': N_imf_err_corrected,
    
    # Integrated masses
    'integrated_mass_corrected': integrated_mass_corrected,
    'integrated_mass_uncorrected': integrated_mass_uncorrected,
    'mean_stellar_mass_corrected': integrated_mass_corrected / total_mass_corrected,
    'mean_stellar_mass_uncorrected': integrated_mass_uncorrected / total_mass_uncorrected,
    
    # Isochrone info
    'isochrone_age_gyr': age_gyr_val,
    'isochrone_metallicity_feh': median_MH,
}

# Save as CSV summary
summary_df = pd.DataFrame({
    'Parameter': [
        'Cluster Name',
        'Age (Gyr)',
        'Distance (pc)',
        'Parallax (mas)',
        'Total Members (uncorrected)',
        'Total Members (corrected)',
        'Completeness Correction Factor',
        'Integrated Mass - Uncorrected (M☉)',
        'Integrated Mass - Corrected (M☉)',
        'Mean Stellar Mass - Uncorrected (M☉)',
        'Mean Stellar Mass - Corrected (M☉)',
    ],
    'Value': [
        cluster_name,
        f'{age_gyr:.1f}',
        f'{distance_pc:.0f}',
        f'{1000.0/distance_pc:.4f}',
        f'{export_data["total_members_uncorrected"]:.0f}',
        f'{export_data["total_members_corrected"]:.0f}',
        f'{correction_factor:.2f}×',
        f'{integrated_mass_uncorrected:.1f}',
        f'{integrated_mass_corrected:.1f}',
        f'{export_data["mean_stellar_mass_uncorrected"]:.3f}',
        f'{export_data["mean_stellar_mass_corrected"]:.3f}',
    ]
})

summary_file = f'Data/{cluster_name.lower()}_mass_summary.csv'
summary_df.to_csv(summary_file, index=False)
print(f"\n✓ Saved summary to: {summary_file}\n{summary_df.to_string(index=False)}")

# Save detailed arrays as numpy
np.savez_compressed(
    f'Data/{cluster_name.lower()}_imf_detailed.npz',
    mass_centers=mass_centers_imf,
    N_imf_corrected=N_imf_corrected,
    N_imf_uncorrected=N_imf_uncorrected,
    N_imf_err_corrected=N_imf_err_corrected,
    r_centers=r_centers,
    sigma_members_corrected=sigma_members_corrected,
    sigma_err_corrected=sigma_err_corrected,
    mag_centers=mag_centers,
    N_lf_corrected=N_lf_corrected,
)
print(f"✓ Saved detailed data to: Data/{cluster_name.lower()}_imf_detailed.npz")

print("\nTo use in mass_estimation.ipynb:")
print("─" * 70)
print("""
# Load summary
import pandas as pd
summary = pd.read_csv(f'Data/{cluster_name}_mass_summary.csv')

# Load detailed arrays
data = np.load(f'Data/{cluster_name}_imf_detailed.npz')
mass_centers = data['mass_centers']
N_imf_corrected = data['N_imf_corrected']
# ... etc
""")
```

---

## What Gets Exported

### Summary Table (CSV)
- Cluster name, age, distance
- Member counts (corrected vs uncorrected)
- Total cluster masses
- Mean stellar masses
- Metallicity

### Detailed Arrays (NPZ)
- Mass centers and IMF counts
- Magnitude centers and luminosity function
- Radial centers and surface density
- Uncertainties on all quantities

---

## Example: Using in mass_estimation.ipynb

Add this cell after importing mass_estimation utilities:

```python
# Load completeness-corrected results from membership.ipynb
CLUSTER = "M2"  # or "M34"

# Load summary
summary = pd.read_csv(f'Data/{cluster_name.lower()}_mass_summary.csv')
cluster_mass_corrected = float(summary[summary['Parameter'] == 'Integrated Mass - Corrected (M☉)']['Value'].values[0])

# Load detailed arrays
data = np.load(f'Data/{CLUSTER.lower()}_imf_detailed.npz')
mass_imf = data['mass_centers']
N_imf_obs = data['N_imf_corrected']

print(f"Cluster: {CLUSTER}")
print(f"Corrected integrated mass: {cluster_mass_corrected:.1f} M☉")
print(f"IMF data shape: {N_imf_obs.shape}")

# Now can fit IMF slope, compute uncertainties, compare to models, etc.
```

---

## Files Created

After running export cell, you'll have:

```
Data/
├── m2_mass_summary.csv                    # Quick reference table
├── m2_imf_detailed.npz                    # All detailed arrays
├── m34_mass_summary.csv
├── m34_imf_detailed.npz
└── completeness_model_params.txt          # (from AST)
```

---

## Advanced: Compare Both Clusters

```python
# Load both M2 and M34
m2_data = np.load('Data/m2_imf_detailed.npz')
m34_data = np.load('Data/m34_imf_detailed.npz')

m2_summary = pd.read_csv('Data/m2_mass_summary.csv')
m34_summary = pd.read_csv('Data/m34_mass_summary.csv')

# Compare
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# IMF comparison
ax = axes[0]
ax.loglog(m2_data['mass_centers'], m2_data['N_imf_corrected'], 
         'o-', label=f"M2 ({m2_summary.iloc[1, 1]})", linewidth=2)
ax.loglog(m34_data['mass_centers'], m34_data['N_imf_corrected'], 
         's-', label=f"M34 ({m34_summary.iloc[1, 1]})", linewidth=2)
ax.set_xlabel('Mass (M☉)')
ax.set_ylabel('Number of Stars')
ax.legend()
ax.grid(alpha=0.3, which='both')

# Radial profile comparison
ax = axes[1]
ax.semilogy(m2_data['r_centers'], m2_data['sigma_members_corrected'], 
           'o-', label='M2', linewidth=2)
ax.semilogy(m34_data['r_centers'], m34_data['sigma_members_corrected'], 
           's-', label='M34', linewidth=2)
ax.set_xlabel('Radius (arcmin)')
ax.set_ylabel('Surface Density (stars/arcmin²)')
ax.legend()
ax.grid(alpha=0.3, which='both')

plt.tight_layout()
plt.show()
```

---

## Notes

1. **File paths** assume `Data/` subdirectory exists
2. **Cluster name** is automatically set based on `CLUSTER` variable in cell 1
3. **Numpy compatibility** – Uses `np.savez_compressed()` for efficient storage
4. **Pandas required** – For CSV export (already imported in membership.ipynb)

---

**Last Updated:** December 8, 2025
