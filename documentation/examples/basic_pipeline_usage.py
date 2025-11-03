"""
Basic Pipeline Usage Example
=============================

This script demonstrates how to run the complete SeaExplorer processing pipeline
from raw data to quality-controlled outputs.

Author: Benedetta Torelli
Repository: https://github.com/BennyTorelli/Seaexplorer_pyglider
Date: October 2025
"""

import sys
sys.path.append('/Users/benedettatorelli/Desktop/PyGlider_SeaExplorer_Project')

from MASTER_pyglider_pipeline import SeaExplorerProcessor
import os

# ============================================================================
# STEP 1: Initialize Processor
# ============================================================================

print("=" * 70)
print("SEAEXPLORER DATA PROCESSING - BASIC EXAMPLE")
print("=" * 70)

# Create processor instance
processor = SeaExplorerProcessor()

# Verify raw data is present
raw_files = len([f for f in os.listdir('input/raw') if f.endswith(('.pld1', '.gli'))])
print(f"\nRaw files detected: {raw_files}")

if raw_files == 0:
    print("ERROR: No raw files found in input/raw/")
    print("Please copy .pld1 and .gli files to input/raw/ directory")
    sys.exit(1)

# ============================================================================
# STEP 2: Run Complete Pipeline
# ============================================================================

print("\nStarting complete pipeline...")
print("This will take approximately 20-30 minutes.\n")

try:
    # Run all processing steps (0 through 6)
    processor.run_complete_pipeline()
    
    print("\n" + "=" * 70)
    print("✓ PIPELINE COMPLETED SUCCESSFULLY!")
    print("=" * 70)
    
except Exception as e:
    print(f"\n✗ Pipeline failed: {e}")
    sys.exit(1)

# ============================================================================
# STEP 3: Run Quality Control
# ============================================================================

print("\nRunning quality control analysis...")

import subprocess

try:
    result = subprocess.run(
        ['python', 'scripts/qc_variables.py', '--latest'],
        capture_output=True,
        text=True,
        check=True
    )
    print(result.stdout)
    print("\n✓ Quality control completed")
    
except subprocess.CalledProcessError as e:
    print(f"✗ QC failed: {e}")
    print(e.stderr)

# ============================================================================
# STEP 4: Verify Outputs
# ============================================================================

print("\n" + "=" * 70)
print("OUTPUT VERIFICATION")
print("=" * 70)

output_dir = 'output/analysis'
expected_files = [
    'sea074-2025_standard_names.nc',
    'seaexplorer_data_complete_standard_names.csv',
    'seaexplorer_qc_variables.csv'
]

for filename in expected_files:
    filepath = os.path.join(output_dir, filename)
    if os.path.exists(filepath):
        size_mb = os.path.getsize(filepath) / (1024 * 1024)
        print(f"✓ {filename:<50} ({size_mb:.1f} MB)")
    else:
        print(f"✗ {filename:<50} (MISSING)")

# ============================================================================
# STEP 5: Basic Data Summary
# ============================================================================

print("\n" + "=" * 70)
print("DATA SUMMARY")
print("=" * 70)

import xarray as xr
import numpy as np

try:
    # Load NetCDF
    ds = xr.open_dataset('output/analysis/sea074-2025_standard_names.nc')
    
    print(f"\nDataset dimensions:")
    print(f"  Total measurements: {len(ds['TIME']):,}")
    print(f"  Time range: {ds['TIME'].values[0]} to {ds['TIME'].values[-1]}")
    
    print(f"\nSpatial coverage:")
    print(f"  Latitude:  {float(ds['LATITUDE'].min()):.4f}° to {float(ds['LATITUDE'].max()):.4f}°N")
    print(f"  Longitude: {float(ds['LONGITUDE'].min()):.4f}° to {float(ds['LONGITUDE'].max()):.4f}°E")
    print(f"  Depth:     {float(ds['DEPTH'].min()):.1f} to {float(ds['DEPTH'].max()):.1f} m")
    
    print(f"\nVariable ranges:")
    print(f"  Temperature:  {float(ds['TEMP'].min()):.2f} to {float(ds['TEMP'].max()):.2f} °C")
    print(f"  Salinity:     {float(ds['PSAL'].min()):.2f} to {float(ds['PSAL'].max()):.2f} PSU")
    print(f"  Oxygen:       {float(ds['DOXY'].min()):.1f} to {float(ds['DOXY'].max()):.1f} µmol/kg")
    print(f"  Chlorophyll:  {float(ds['CHLA'].min()):.3f} to {float(ds['CHLA'].max()):.3f} mg/m³")
    
    ds.close()
    
except Exception as e:
    print(f"Could not load dataset: {e}")

print("\n" + "=" * 70)
print("PROCESSING COMPLETE!")
print("=" * 70)
print("\nNext steps:")
print("  1. Review QC results in: output/analysis/seaexplorer_qc_variables.csv")
print("  2. Visualize data with QGIS or ODV")
print("  3. Run custom analyses on NetCDF or CSV outputs")
print("\nDocumentation: https://github.com/BennyTorelli/Seaexplorer_pyglider/tree/main/documentation")
