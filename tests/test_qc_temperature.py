import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import xarray as xr
from scripts.qc_variables import range_qc_temperature, export_qc_to_csv


def test_range_qc_all_pass():
    time = np.arange(3).astype('datetime64[ns]')
    # Use values within LaPalma range [5, 30]°C
    temp = xr.DataArray([10.0, 15.0, 25.0], coords=[time], dims=['time'])
    ds = xr.Dataset({'temperature': temp})
    ds2 = range_qc_temperature(ds)
    assert 'temperature_qc' in ds2
    assert np.array_equal(ds2['temperature_qc'].values, np.array([1, 1, 1], dtype=np.int8))


def test_range_qc_fail_and_missing():
    time = np.arange(4).astype('datetime64[ns]')
    # LaPalma range [5, 30]°C: one too low (0), one nan, one too high (35), one ok (20)
    tempvals = np.array([0.0, np.nan, 35.0, 20.0], dtype=float)
    temp = xr.DataArray(tempvals, coords=[time], dims=['time'])
    ds = xr.Dataset({'temperature': temp})
    ds2 = range_qc_temperature(ds)
    flags = ds2['temperature_qc'].values
    assert flags[0] == 4  # too low
    assert np.isnan(flags[1])  # missing - no flag (NaN)
    assert flags[2] == 4  # too high
    assert flags[3] == 1  # ok


def test_export_csv(tmp_path):
    time = np.arange(3).astype('datetime64[ns]')
    # LaPalma range [5, 30]°C: 10 (pass), 3 (fail - too low), nan (missing - no flag)
    tempvals = np.array([10.0, 3.0, np.nan], dtype=float)
    temp = xr.DataArray(tempvals, coords=[time], dims=['time'])
    
    # Add LATITUDE and LONGITUDE (realistic glider coordinates)
    lat = xr.DataArray([28.5, 28.6, 28.7], coords=[time], dims=['time'])  # Near La Palma
    lon = xr.DataArray([-17.8, -17.9, -18.0], coords=[time], dims=['time'])  # Near La Palma
    
    ds = xr.Dataset({'temperature': temp, 'LATITUDE': lat, 'LONGITUDE': lon})
    ds2 = range_qc_temperature(ds)
    csv_file = tmp_path / 'out.csv'
    out = export_qc_to_csv(ds2, csv_path=str(csv_file))
    assert out == str(csv_file)
    import pandas as pd
    df = pd.read_csv(out)
    
    # Check Range_QC column
    assert 'TEMP_Range_QC' in df.columns
    assert df['TEMP_Range_QC'].iloc[0] == 1  # pass
    assert df['TEMP_Range_QC'].iloc[1] == 4  # fail
    assert pd.isna(df['TEMP_Range_QC'].iloc[2])  # missing - no range flag
    
    # Check Na_QC column (missing value flags)
    assert 'TEMP_Na_QC' in df.columns
    assert df['TEMP_Na_QC'].iloc[0] == 1  # value present
    assert df['TEMP_Na_QC'].iloc[1] == 1  # value present
    assert df['TEMP_Na_QC'].iloc[2] == 9  # value missing
    
    # Check that Location_QC column exists (combined LATITUDE and LONGITUDE QC)
    assert 'Location_QC' in df.columns
    # All coordinates should pass (valid range)
    assert all(df['Location_QC'] == 1)


