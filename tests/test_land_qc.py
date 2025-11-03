#!/usr/bin/env python3
"""
Test script for LAND_QC functionality.
Creates a test NetCDF file with mixed coordinates (some on land, some in ocean)
and verifies that the QC script correctly identifies land points.
"""

import xarray as xr
import numpy as np
import pandas as pd
import os
import sys

# Add parent directory to path to import qc_variables
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def create_test_netcdf():
    """
    Create a test NetCDF file with various coordinate scenarios:
    1. Points clearly in the ocean (should get flag 1 - GOOD)
    2. Points on La Palma island (should get flag 4 - BAD)
    3. Points near the coast but in ocean (should get flag 1 - GOOD)
    """
    
    # Define test coordinates
    # La Palma bounds: Lon -18.01 to -17.72, Lat 28.45 to 28.86
    
    test_points = [
        # Description, Latitude, Longitude, Expected_Flag
        ("Ocean - West of La Palma", 28.65, -18.10, 1),
        ("Ocean - East of La Palma", 28.65, -17.70, 1),
        ("Ocean - South of La Palma", 28.40, -17.85, 1),
        ("Ocean - North of La Palma", 28.90, -17.85, 1),
        ("LAND - Center of La Palma", 28.65, -17.85, 4),
        ("LAND - Santa Cruz de La Palma", 28.68, -17.76, 4),
        ("LAND - Los Llanos", 28.66, -17.91, 4),
        ("LAND - Fuencaliente (south tip)", 28.50, -17.83, 4),
        ("LAND - Punta Cumplida (north)", 28.82, -17.75, 4),
        ("Ocean - Near coast but outside", 28.70, -18.02, 1),
    ]
    
    n_points = len(test_points)
    
    # Create time array
    time = pd.date_range('2025-02-12', periods=n_points, freq='1H')
    
    # Extract coordinates and descriptions
    descriptions = [p[0] for p in test_points]
    lats = np.array([p[1] for p in test_points])
    lons = np.array([p[2] for p in test_points])
    expected_flags = np.array([p[3] for p in test_points])
    
    # Create dummy data for other variables (to match real NetCDF structure)
    depth = np.linspace(0, 100, n_points)
    temp = np.full(n_points, 20.0)
    cndc = np.full(n_points, 4.5)
    pres = np.linspace(0, 100, n_points)
    chla = np.full(n_points, 0.5)
    turb = np.full(n_points, 0.1)
    doxy = np.full(n_points, 220.0)
    psal = np.full(n_points, 36.0)
    
    # Create xarray Dataset
    ds = xr.Dataset(
        {
            'TEMP': (['time'], temp, {'units': 'degrees_C', 'long_name': 'Temperature'}),
            'CNDC': (['time'], cndc, {'units': 'S/m', 'long_name': 'Conductivity'}),
            'PRES': (['time'], pres, {'units': 'dbar', 'long_name': 'Pressure'}),
            'CHLA': (['time'], chla, {'units': 'mg/m^3', 'long_name': 'Chlorophyll-a'}),
            'TURB': (['time'], turb, {'units': 'NTU', 'long_name': 'Turbidity'}),
            'DOXY': (['time'], doxy, {'units': 'umol/L', 'long_name': 'Dissolved Oxygen'}),
            'PSAL': (['time'], psal, {'units': 'PSU', 'long_name': 'Practical Salinity'}),
            'LATITUDE': (['time'], lats, {'units': 'degrees_north', 'long_name': 'Latitude'}),
            'LONGITUDE': (['time'], lons, {'units': 'degrees_east', 'long_name': 'Longitude'}),
            'depth': (['time'], depth, {'units': 'm', 'long_name': 'Depth'}),
        },
        coords={
            'time': time,
        },
        attrs={
            'title': 'Test NetCDF for LAND_QC validation',
            'institution': 'PyGlider Test Suite',
            'source': 'Synthetic test data',
        }
    )
    
    # Save to file
    output_dir = 'tests/data/'
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, 'test_land_qc.nc')
    ds.to_netcdf(output_file)
    
    print(f"✓ Created test NetCDF: {output_file}")
    print(f"\n{'='*80}")
    print(f"TEST COORDINATES:")
    print(f"{'='*80}")
    for i, desc in enumerate(descriptions):
        print(f"{i+1}. {desc}")
        print(f"   LAT: {lats[i]:.6f}, LON: {lons[i]:.6f}")
        print(f"   Expected LAND_QC: {expected_flags[i]} ({'GOOD - Ocean' if expected_flags[i] == 1 else 'BAD - Land'})")
    print(f"{'='*80}\n")
    
    # Save expected results for validation
    expected_df = pd.DataFrame({
        'point': range(1, n_points + 1),
        'description': descriptions,
        'latitude': lats,
        'longitude': lons,
        'expected_land_qc': expected_flags
    })
    expected_file = os.path.join(output_dir, 'expected_land_qc.csv')
    expected_df.to_csv(expected_file, index=False)
    print(f"✓ Saved expected results: {expected_file}\n")
    
    return output_file, expected_df


def run_qc_script(test_nc_file):
    """Run the qc_variables.py script on the test NetCDF file."""
    from scripts.qc_variables import export_qc_to_csv
    
    # Open test NetCDF
    ds = xr.open_dataset(test_nc_file)
    
    # Apply range QC (required before export)
    from scripts.qc_variables import range_qc_temperature
    ds = range_qc_temperature(ds, varname='TEMP')
    
    # Export to CSV (this will compute LAND_QC)
    output_csv = 'tests/data/test_land_qc_results.csv'
    export_qc_to_csv(ds, varname='TEMP', csv_path=output_csv)
    
    print(f"✓ QC analysis completed: {output_csv}\n")
    
    ds.close()
    return output_csv


def validate_results(results_csv, expected_df):
    """Compare actual LAND_QC results with expected values."""
    
    # Read results
    results = pd.read_csv(results_csv)
    
    print(f"{'='*80}")
    print(f"VALIDATION RESULTS:")
    print(f"{'='*80}\n")
    
    # Extract LAND_QC column
    land_qc = results['LAND_QC'].values
    expected = expected_df['expected_land_qc'].values
    
    # Compare
    all_correct = True
    for i in range(len(expected)):
        actual = land_qc[i]
        exp = expected[i]
        status = "✓ PASS" if actual == exp else "✗ FAIL"
        
        desc = expected_df['description'].iloc[i]
        lat = expected_df['latitude'].iloc[i]
        lon = expected_df['longitude'].iloc[i]
        
        print(f"{i+1}. {desc}")
        print(f"   LAT: {lat:.6f}, LON: {lon:.6f}")
        print(f"   Expected: {exp} | Actual: {actual} | {status}")
        
        if actual != exp:
            all_correct = False
            print(f"   ⚠️  ERROR: Flag mismatch!")
        print()
    
    print(f"{'='*80}")
    if all_correct:
        print(f"✓✓✓ ALL TESTS PASSED! ✓✓✓")
        print(f"LAND_QC is working correctly!")
    else:
        print(f"✗✗✗ SOME TESTS FAILED ✗✗✗")
        print(f"Please check the implementation.")
    print(f"{'='*80}\n")
    
    return all_correct


if __name__ == '__main__':
    print("\n" + "="*80)
    print("LAND_QC TEST SUITE")
    print("="*80 + "\n")
    
    # Step 1: Create test NetCDF
    print("STEP 1: Creating test NetCDF with mixed coordinates...")
    test_nc, expected_df = create_test_netcdf()
    
    # Step 2: Run QC script
    print("STEP 2: Running qc_variables.py on test data...")
    results_csv = run_qc_script(test_nc)
    
    # Step 3: Validate results
    print("STEP 3: Validating LAND_QC results...")
    success = validate_results(results_csv, expected_df)
    
    # Exit with appropriate code
    sys.exit(0 if success else 1)
