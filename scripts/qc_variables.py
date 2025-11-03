"""QC utilities for multiple variables (TEMP, CNDC, ...)."""

from typing import Tuple, Optional
import xarray as xr
import numpy as np
import pandas as pd
import logging
import json
from shapely.geometry import Point, Polygon

_log = logging.getLogger(__name__)

# QC Range constants for La Palma zone (SeaExplorer sensors)
# Source: Coriolis/PLOCAN QC ranges for oceanographic data
QC_RANGES_LAPALMA = {
    'TEMP': {'min': 5.0, 'max': 30.0, 'unit': '°C'},
    'PSAL': {'min': 33.0, 'max': 38.0, 'unit': 'PSU'},
    'CNDC': {'min': 0.0, 'max': 7.0, 'unit': 'S/m'},  # Derived from PSAL range
    'DOXY': {'min': 110.0, 'max': 250.0, 'unit': 'µmol/L'},
    'CHLA': {'min': 0.0, 'max': 42.0, 'unit': 'mg/m³'},
    'TURB': {'min': 0.0, 'max': 5.0, 'unit': 'NTU'},
    'PRES': {'min': 0.0, 'max': 2000.0, 'unit': 'dbar'},  # Typical glider max depth
}


def range_qc_temperature(ds: xr.Dataset, varname: Optional[str] = None,
                         valid_min: float = None, valid_max: float = None) -> xr.Dataset:
    """Apply range QC to a temperature variable in an xarray.Dataset.

        Adds or updates a companion quality flag variable named `<varname>_qc` with
        integer flags following the requested convention:
            1 = PASS (in-range)
            4 = FAIL (out-of-range)
            9 = MISSING (NaN)
            0 = NOT EVALUATED

    The function preserves attributes on the original variable and adds
    reasonable attributes to the QC variable.

    Args:
        ds: xarray.Dataset containing the variable to check.
        varname: variable name to QC (default: 'temperature').
        valid_min: minimum valid value (inclusive). If None, uses QC_RANGES_LAPALMA['TEMP']['min'].
        valid_max: maximum valid value (inclusive). If None, uses QC_RANGES_LAPALMA['TEMP']['max'].

    Returns:
        The dataset with an added `<varname>_qc` variable (int8).
    """
    # Use LaPalma ranges as defaults
    if valid_min is None:
        valid_min = QC_RANGES_LAPALMA['TEMP']['min']
    if valid_max is None:
        valid_max = QC_RANGES_LAPALMA['TEMP']['max']
    
    # Backwards-compatible: accept either standardized 'TEMP' or legacy 'temperature'
    if varname is None:
        if 'TEMP' in ds:
            varname = 'TEMP'
            _log.info("Using standardized temperature variable '%s'", varname)
        elif 'temperature' in ds:
            varname = 'temperature'
            _log.info("Using legacy temperature variable '%s'", varname)
        else:
            raise KeyError("Variable 'TEMP' or 'temperature' not found in dataset; this QC operates only on temperature variables")
    else:
        # allow either accepted name
        if varname not in ('TEMP', 'temperature'):
            raise KeyError("This QC function only accepts the standardized variable name 'TEMP' or legacy 'temperature'")

    data = ds[varname].values
    # initialize flags as NaN (no flag for missing values)
    flags = np.full(data.shape, np.nan, dtype=float)

    # Where values are finite, mark PASS or FAIL
    finite_mask = np.isfinite(data)
    if finite_mask.any():
        vals = data[finite_mask]
        # Only evaluate present values: pass=1, fail=4
        pass_mask = (vals >= valid_min) & (vals <= valid_max)
        flags[finite_mask] = np.where(pass_mask, 1, 4)

    qc_name = f"{varname}_qc"
    attrs = {
        'long_name': f'{varname} quality flag',
        'flag_values': np.array([1, 4], dtype=np.int8),
        'flag_meanings': 'pass fail',
        'units': '1',
        'method': 'range check',
        'valid_min': 1,
        'valid_max': 4,
        'comment': f'Valid range: {valid_min} to {valid_max}. Flag 4 indicates out of range. Missing values have no flag.',
    }

    # Insert or replace qc variable
    # preserve dataset variable name casing
    ds[qc_name] = (ds[varname].dims, flags, attrs)
    _log.info('Applied range QC to %s (%g..%g)', varname, valid_min, valid_max)
    return ds


def range_qc_variable(ds: xr.Dataset, varname: str,
                      valid_min: float, valid_max: float) -> xr.Dataset:
    """Generic range QC for any existing variable name (expects standardized name).

    Adds `<varname>_qc` to the dataset with flags 0=not evaluated, 1=pass, 4=fail.
    Missing values (NaN) get flag 0 (not evaluated).
    """
    if varname not in ds:
        raise KeyError(f"Variable '{varname}' not found in dataset")

    data = ds[varname].values
    flags = np.zeros(data.shape, dtype=float)  # 0 for not evaluated (including missing)
    finite_mask = np.isfinite(data)
    if finite_mask.any():
        vals = data[finite_mask]
        # Only evaluate present values: pass=1, fail=4
        pass_mask = (vals >= valid_min) & (vals <= valid_max)
        flags[finite_mask] = np.where(pass_mask, 1, 4)

    qc_name = f"{varname}_qc"
    attrs = {
        'long_name': f'{varname} quality flag',
        'flag_values': np.array([0, 1, 4], dtype=np.int8),
        'flag_meanings': 'not_evaluated pass fail',
        'units': '1',
        'method': 'range check',
        'valid_min': 0,
        'valid_max': 4,
        'comment': f'Valid range: {valid_min} to {valid_max}. Flag 0 = not evaluated (missing data), Flag 1 = pass, Flag 4 = out of range.',
    }
    ds[qc_name] = (ds[varname].dims, flags, attrs)
    _log.info('Applied range QC to %s (%g..%g)', varname, valid_min, valid_max)
    return ds


def range_qc_time(ds: xr.Dataset, time_coord: str = 'TIME') -> xr.Dataset:
    """Apply range QC to time coordinate.

    Checks if timestamps are within valid range: 1950-01-01 to current date.
    
    Args:
        ds: xarray.Dataset containing the time coordinate.
        time_coord: name of the time coordinate (default: 'TIME').
    
    Returns:
        The dataset with an added `TIME_qc` variable (int8).
    """
    import pandas as pd
    from datetime import datetime
    
    # Find time coordinate (case-insensitive)
    time_var = None
    for coord in ds.coords:
        if coord.upper() == time_coord.upper():
            time_var = coord
            break
    
    if time_var is None:
        raise KeyError(f"Time coordinate '{time_coord}' not found in dataset")
    
    # Define valid time range
    min_time = pd.Timestamp('1950-01-01')
    max_time = pd.Timestamp(datetime.now())  # Current date at execution time
    
    time_values = pd.to_datetime(ds[time_var].values)
    flags = np.full(time_values.shape, np.nan, dtype=float)  # NaN for missing
    
    # Only evaluate non-missing timestamps
    valid_mask = pd.notna(time_values)
    if valid_mask.any():
        valid_times = time_values[valid_mask]
        # Check range: pass=1, fail=4 (only for present values)
        pass_mask = (valid_times >= min_time) & (valid_times <= max_time)
        flags[valid_mask] = np.where(pass_mask, 1, 4)
    
    qc_name = 'TIME_qc'
    attrs = {
        'long_name': 'time quality flag',
        'flag_values': np.array([1, 4], dtype=np.int8),
        'flag_meanings': 'pass fail',
        'units': '1',
        'method': 'range check',
        'valid_min': 1,
        'valid_max': 4,
        'comment': f'Valid range: 1950-01-01 to {max_time.strftime("%Y-%m-%d")}. Flag 4 indicates out of range. Missing values have no flag.',
    }
    
    ds[qc_name] = (ds[time_var].dims, flags, attrs)
    _log.info('Applied time QC: valid range 1950-01-01 to %s', max_time.strftime('%Y-%m-%d'))
    return ds


def range_qc_latitude(ds: xr.Dataset, lat_var: str = 'LATITUDE') -> xr.Dataset:
    """Apply range QC to latitude coordinate.

    Checks if latitude values are within valid range: -90 to 90 degrees.
    
    Args:
        ds: xarray.Dataset containing the latitude variable.
        lat_var: name of the latitude variable (default: 'LATITUDE').
    
    Returns:
        The dataset with an added `LATITUDE_qc` variable (int8).
    """
    # Find latitude variable (case-insensitive)
    latitude_var = None
    for var in list(ds.data_vars) + list(ds.coords):
        if var.upper() == lat_var.upper():
            latitude_var = var
            break
    
    if latitude_var is None:
        raise KeyError(f"Latitude variable '{lat_var}' not found in dataset")
    
    # Define valid latitude range
    min_lat = -90.0
    max_lat = 90.0
    
    lat_values = ds[latitude_var].values
    flags = np.full(lat_values.shape, np.nan, dtype=float)  # NaN for missing
    
    # Only evaluate non-missing values
    valid_mask = ~np.isnan(lat_values)
    if valid_mask.any():
        valid_lats = lat_values[valid_mask]
        # Check range: pass=1, fail=4 (only for present values)
        pass_mask = (valid_lats >= min_lat) & (valid_lats <= max_lat)
        flags[valid_mask] = np.where(pass_mask, 1, 4)
    
    qc_name = 'LATITUDE_qc'
    attrs = {
        'long_name': 'latitude quality flag',
        'flag_values': np.array([1, 4], dtype=np.int8),
        'flag_meanings': 'pass fail',
        'units': '1',
        'method': 'range check',
        'valid_min': 1,
        'valid_max': 4,
        'comment': f'Valid range: {min_lat} to {max_lat} degrees. Flag 4 indicates out of range. Missing values have no flag.',
    }
    
    ds[qc_name] = (ds[latitude_var].dims, flags, attrs)
    _log.info('Applied latitude QC: valid range %.1f to %.1f degrees', min_lat, max_lat)
    return ds


def range_qc_longitude(ds: xr.Dataset, lon_var: str = 'LONGITUDE') -> xr.Dataset:
    """Apply range QC to longitude coordinate.

    Checks if longitude values are within valid range: -180 to 180 degrees.
    
    Args:
        ds: xarray.Dataset containing the longitude variable.
        lon_var: name of the longitude variable (default: 'LONGITUDE').
    
    Returns:
        The dataset with an added `LONGITUDE_qc` variable (int8).
    """
    # Find longitude variable (case-insensitive)
    longitude_var = None
    for var in list(ds.data_vars) + list(ds.coords):
        if var.upper() == lon_var.upper():
            longitude_var = var
            break
    
    if longitude_var is None:
        raise KeyError(f"Longitude variable '{lon_var}' not found in dataset")
    
    # Define valid longitude range
    min_lon = -180.0
    max_lon = 180.0
    
    lon_values = ds[longitude_var].values
    flags = np.full(lon_values.shape, np.nan, dtype=float)  # NaN for missing
    
    # Only evaluate non-missing values
    valid_mask = ~np.isnan(lon_values)
    if valid_mask.any():
        valid_lons = lon_values[valid_mask]
        # Check range: pass=1, fail=4 (only for present values)
        pass_mask = (valid_lons >= min_lon) & (valid_lons <= max_lon)
        flags[valid_mask] = np.where(pass_mask, 1, 4)
    
    qc_name = 'LONGITUDE_qc'
    attrs = {
        'long_name': 'longitude quality flag',
        'flag_values': np.array([1, 4], dtype=np.int8),
        'flag_meanings': 'pass fail',
        'units': '1',
        'method': 'range check',
        'valid_min': 1,
        'valid_max': 4,
        'comment': f'Valid range: {min_lon} to {max_lon} degrees. Flag 4 indicates out of range. Missing values have no flag.',
    }
    
    ds[qc_name] = (ds[longitude_var].dims, flags, attrs)
    _log.info('Applied longitude QC: valid range %.1f to %.1f degrees', min_lon, max_lon)
    return ds


def export_qc_to_csv(ds: xr.Dataset, varname: str = None, csv_path: str = None):
    """Export a compact QC CSV containing only time, depth, and QC/temperature.

    The output columns will be exactly: time, depth, QC/temperature
    where QC/temperature uses the lowercase legacy label you requested.
    The function auto-detects `TEMP` vs `temperature` when `varname` is None.
    """
    # Apply TIME QC first
    ds = range_qc_time(ds)
    
    # Apply LATITUDE and LONGITUDE QC if variables exist
    try:
        ds = range_qc_latitude(ds)
    except KeyError:
        _log.info("LATITUDE not found, skipping latitude QC")
    
    try:
        ds = range_qc_longitude(ds)
    except KeyError:
        _log.info("LONGITUDE not found, skipping longitude QC")
    
    # Accept standardized 'TEMP' or legacy 'temperature'
    if varname is None:
        if 'TEMP' in ds:
            varname = 'TEMP'
            _log.info("Using standardized temperature variable '%s' for CSV export", varname)
        elif 'temperature' in ds:
            varname = 'temperature'
            _log.info("Using legacy temperature variable '%s' for CSV export", varname)
        else:
            raise KeyError("Variable 'TEMP' or 'temperature' not found in dataset; export requires a temperature variable")
    qc_name = f"{varname}_qc"
    if qc_name not in ds:
        raise KeyError(f"QC variable '{qc_name}' not found in dataset. Run range_qc_temperature first")

    # Build a dataframe from the qc variable and bring in the depth coordinate if present
    df_qc = ds[qc_name].to_dataframe(name='qc').reset_index()

    # Normalize column names: time may be 'TIME' in standardized files
    # find the time column present in the dataframe
    time_col = None
    for c in df_qc.columns:
        if c.lower() == 'time':
            time_col = c
            break
    if time_col is None:
        # fallback: use the first column as time-like
        time_col = df_qc.columns[0]

    # depth column: try common names in dataset coords
    depth_vals = None
    for coord in ('DEPTH', 'depth'):
        if coord in ds.coords:
            depth_vals = ds[coord].values
            break

    # Build compact dataframe for TEMP
    compact = df_qc[[time_col, 'qc']].copy()
    compact = compact.rename(columns={time_col: 'time', 'qc': 'TEMP_Range_QC'})
    # attach depth if available (align by index length)
    if depth_vals is not None:
        # depth might be a per-sample array aligned with the dataset
        compact['depth'] = depth_vals
    else:
        compact['depth'] = ''

    # Add TIME QC column (renamed to Date_QC)
    if 'TIME_qc' in ds:
        time_qc = ds['TIME_qc'].to_dataframe(name='time_qc').reset_index()
        if time_col in time_qc.columns:
            time_qc = time_qc[[time_col, 'time_qc']]
            time_qc = time_qc.rename(columns={time_col: 'time', 'time_qc': 'Date_QC'})
            compact = compact.merge(time_qc, on='time', how='left')
        else:
            compact['Date_QC'] = ''
    else:
        compact['Date_QC'] = ''

    # Add LATITUDE and LONGITUDE QC columns (combined as Location_QC)
    # We'll use the LATITUDE_qc for Location_QC (assuming both coordinates are checked together)
    if 'LATITUDE_qc' in ds and 'LONGITUDE_qc' in ds:
        # Combine lat/lon QC: fail if either fails
        lat_qc_vals = ds['LATITUDE_qc'].values
        lon_qc_vals = ds['LONGITUDE_qc'].values
        # Location QC: 1 if both pass, 4 if either fails, NaN if either is NaN
        location_qc = np.full(lat_qc_vals.shape, np.nan, dtype=float)
        both_present = np.isfinite(lat_qc_vals) & np.isfinite(lon_qc_vals)
        if both_present.any():
            # Where both have flags: 1 if both pass, 4 if either fails
            location_qc[both_present] = np.where((lat_qc_vals[both_present] == 4) | (lon_qc_vals[both_present] == 4), 4, 1)
        
        location_df = ds['LATITUDE_qc'].to_dataframe(name='loc_qc').reset_index()
        if time_col in location_df.columns:
            location_df = location_df[[time_col]]
            location_df['Location_QC'] = location_qc
            location_df = location_df.rename(columns={time_col: 'time'})
            compact = compact.merge(location_df, on='time', how='left')
        else:
            compact['Location_QC'] = ''
    elif 'LATITUDE_qc' in ds:
        # Only latitude available
        lat_qc = ds['LATITUDE_qc'].to_dataframe(name='lat_qc').reset_index()
        if time_col in lat_qc.columns:
            lat_qc = lat_qc[[time_col, 'lat_qc']]
            lat_qc = lat_qc.rename(columns={time_col: 'time', 'lat_qc': 'Location_QC'})
            compact = compact.merge(lat_qc, on='time', how='left')
        else:
            compact['Location_QC'] = ''
    elif 'LONGITUDE_qc' in ds:
        # Only longitude available
        lon_qc = ds['LONGITUDE_qc'].to_dataframe(name='lon_qc').reset_index()
        if time_col in lon_qc.columns:
            lon_qc = lon_qc[[time_col, 'lon_qc']]
            lon_qc = lon_qc.rename(columns={time_col: 'time', 'lon_qc': 'Location_QC'})
            compact = compact.merge(lon_qc, on='time', how='left')
        else:
            compact['Location_QC'] = ''
    else:
        compact['Location_QC'] = ''

    # Compute QC for additional variables if present (CNDC, DOXY, CHLA)
    
    # TEMP Sensor QC (physical sensor range -5 to 42°C)
    temp_var = 'TEMP' if 'TEMP' in ds else ('temperature' if 'temperature' in ds else None)
    if temp_var:
        temp_sensor_min = -5.0
        temp_sensor_max = 42.0
        ds_temp_sensor = range_qc_variable(ds, temp_var, temp_sensor_min, temp_sensor_max)
        temp_sensor_q = ds_temp_sensor[f'{temp_var}_qc'].to_dataframe(name='temp_sensor_q').reset_index()
        if time_col in temp_sensor_q.columns:
            temp_sensor_q = temp_sensor_q[[time_col, 'temp_sensor_q']]
            temp_sensor_q = temp_sensor_q.rename(columns={time_col: 'time', 'temp_sensor_q': 'TEMP_Sensor_QC'})
            compact = compact.merge(temp_sensor_q, on='time', how='left')
        else:
            compact['TEMP_Sensor_QC'] = ''
    else:
        compact['TEMP_Sensor_QC'] = ''
    
    # CNDC QC
    if 'CNDC' in ds:
        cndc_min = QC_RANGES_LAPALMA['CNDC']['min']
        cndc_max = QC_RANGES_LAPALMA['CNDC']['max']
        ds = range_qc_variable(ds, 'CNDC', cndc_min, cndc_max)
        cndc_q = ds['CNDC_qc'].to_dataframe(name='cndc_q').reset_index()
        if time_col in cndc_q.columns:
            cndc_q = cndc_q[[time_col, 'cndc_q']]
            cndc_q = cndc_q.rename(columns={time_col: 'time', 'cndc_q': 'CNDC_Range_QC'})
            compact = compact.merge(cndc_q, on='time', how='left')
        else:
            compact['CNDC_Range_QC'] = ''
    else:
        compact['CNDC_Range_QC'] = ''

    # CNDC Sensor QC - check sensor range [0, 8.5] S/m
    if 'CNDC' in ds:
        # Create sensor QC with different range than LaPalma QC
        cndc_sensor_min = 0.0
        cndc_sensor_max = 8.5
        ds_sensor = range_qc_variable(ds, 'CNDC', cndc_sensor_min, cndc_sensor_max)
        cndc_sensor_q = ds_sensor['CNDC_qc'].to_dataframe(name='cndc_sensor_q').reset_index()
        if time_col in cndc_sensor_q.columns:
            cndc_sensor_q = cndc_sensor_q[[time_col, 'cndc_sensor_q']]
            cndc_sensor_q = cndc_sensor_q.rename(columns={time_col: 'time', 'cndc_sensor_q': 'CNDC_Sensor_QC'})
            compact = compact.merge(cndc_sensor_q, on='time', how='left')
        else:
            compact['CNDC_Sensor_QC'] = ''
    else:
        compact['CNDC_Sensor_QC'] = ''

    # DOXY QC
    if 'DOXY' in ds:
        doxy_min = QC_RANGES_LAPALMA['DOXY']['min']
        doxy_max = QC_RANGES_LAPALMA['DOXY']['max']
        ds = range_qc_variable(ds, 'DOXY', doxy_min, doxy_max)
        doxy_q = ds['DOXY_qc'].to_dataframe(name='doxy_q').reset_index()
        if time_col in doxy_q.columns:
            doxy_q = doxy_q[[time_col, 'doxy_q']]
            doxy_q = doxy_q.rename(columns={time_col: 'time', 'doxy_q': 'DOXY_Range_QC'})
            compact = compact.merge(doxy_q, on='time', how='left')
        else:
            compact['DOXY_Range_QC'] = ''
    else:
        compact['DOXY_Range_QC'] = ''

    # DOXY Sensor QC (physical sensor range 0-500 µmol/L)
    if 'DOXY' in ds:
        doxy_sensor_min = 0.0
        doxy_sensor_max = 500.0
        ds_doxy_sensor = range_qc_variable(ds, 'DOXY', doxy_sensor_min, doxy_sensor_max)
        doxy_sensor_q = ds_doxy_sensor['DOXY_qc'].to_dataframe(name='doxy_sensor_q').reset_index()
        if time_col in doxy_sensor_q.columns:
            doxy_sensor_q = doxy_sensor_q[[time_col, 'doxy_sensor_q']]
            doxy_sensor_q = doxy_sensor_q.rename(columns={time_col: 'time', 'doxy_sensor_q': 'DOXY_Sensor_QC'})
            compact = compact.merge(doxy_sensor_q, on='time', how='left')
        else:
            compact['DOXY_Sensor_QC'] = ''
    else:
        compact['DOXY_Sensor_QC'] = ''

    # PSAL QC - check for both 'PSAL' (standardized) and 'salinity' (legacy)
    psal_var = None
    if 'PSAL' in ds:
        psal_var = 'PSAL'
    elif 'salinity' in ds:
        psal_var = 'salinity'
        # Convert salinity from 1e-3 to PSU before QC
        # Create a copy in the dataset with converted values for QC calculation
        ds['salinity_psu'] = ds['salinity'] * 1000
        ds['salinity_psu'].attrs = ds['salinity'].attrs.copy()
        ds['salinity_psu'].attrs['units'] = 'PSU'
        psal_var = 'salinity_psu'
        _log.info("Converted salinity to PSU for QC calculation")
    
    if psal_var:
        psal_min = QC_RANGES_LAPALMA['PSAL']['min']
        psal_max = QC_RANGES_LAPALMA['PSAL']['max']
        ds = range_qc_variable(ds, psal_var, psal_min, psal_max)
        psal_qc_name = f'{psal_var}_qc'
        psal_q = ds[psal_qc_name].to_dataframe(name='psal_q').reset_index()
        if time_col in psal_q.columns:
            psal_q = psal_q[[time_col, 'psal_q']]
            psal_q = psal_q.rename(columns={time_col: 'time', 'psal_q': 'PSAL_Range_QC'})
            compact = compact.merge(psal_q, on='time', how='left')
        else:
            compact['PSAL_Range_QC'] = ''
    else:
        compact['PSAL_Range_QC'] = ''

    # CHLA QC
    if 'CHLA' in ds:
        chla_min = QC_RANGES_LAPALMA['CHLA']['min']
        chla_max = QC_RANGES_LAPALMA['CHLA']['max']
        ds = range_qc_variable(ds, 'CHLA', chla_min, chla_max)
        chla_q = ds['CHLA_qc'].to_dataframe(name='chla_q').reset_index()
        if time_col in chla_q.columns:
            chla_q = chla_q[[time_col, 'chla_q']]
            chla_q = chla_q.rename(columns={time_col: 'time', 'chla_q': 'CHLA_Range_QC'})
            compact = compact.merge(chla_q, on='time', how='left')
        else:
            compact['CHLA_Range_QC'] = ''
    else:
        compact['CHLA_Range_QC'] = ''

    # CHLA Sensor QC (physical sensor range 0-50 mg/m³)
    if 'CHLA' in ds:
        chla_sensor_min = 0.0
        chla_sensor_max = 50.0
        ds_chla_sensor = range_qc_variable(ds, 'CHLA', chla_sensor_min, chla_sensor_max)
        chla_sensor_q = ds_chla_sensor['CHLA_qc'].to_dataframe(name='chla_sensor_q').reset_index()
        if time_col in chla_sensor_q.columns:
            chla_sensor_q = chla_sensor_q[[time_col, 'chla_sensor_q']]
            chla_sensor_q = chla_sensor_q.rename(columns={time_col: 'time', 'chla_sensor_q': 'CHLA_Sensor_QC'})
            compact = compact.merge(chla_sensor_q, on='time', how='left')
        else:
            compact['CHLA_Sensor_QC'] = ''
    else:
        compact['CHLA_Sensor_QC'] = ''

    # TURB QC
    if 'TURB' in ds:
        turb_min = QC_RANGES_LAPALMA['TURB']['min']
        turb_max = QC_RANGES_LAPALMA['TURB']['max']
        ds = range_qc_variable(ds, 'TURB', turb_min, turb_max)
        turb_q = ds['TURB_qc'].to_dataframe(name='turb_q').reset_index()
        if time_col in turb_q.columns:
            turb_q = turb_q[[time_col, 'turb_q']]
            turb_q = turb_q.rename(columns={time_col: 'time', 'turb_q': 'TURB_Range_QC'})
            compact = compact.merge(turb_q, on='time', how='left')
        else:
            compact['TURB_Range_QC'] = ''
    else:
        compact['TURB_Range_QC'] = ''

    # TURB Sensor QC (physical sensor range 0-25 NTU)
    if 'TURB' in ds:
        turb_sensor_min = 0.0
        turb_sensor_max = 25.0
        ds_turb_sensor = range_qc_variable(ds, 'TURB', turb_sensor_min, turb_sensor_max)
        turb_sensor_q = ds_turb_sensor['TURB_qc'].to_dataframe(name='turb_sensor_q').reset_index()
        if time_col in turb_sensor_q.columns:
            turb_sensor_q = turb_sensor_q[[time_col, 'turb_sensor_q']]
            turb_sensor_q = turb_sensor_q.rename(columns={time_col: 'time', 'turb_sensor_q': 'TURB_Sensor_QC'})
            compact = compact.merge(turb_sensor_q, on='time', how='left')
        else:
            compact['TURB_Sensor_QC'] = ''
    else:
        compact['TURB_Sensor_QC'] = ''

    # LAND_QC: Check if glider coordinates are inside La Palma island polygon
    # Flag 4 (BAD) = inside island polygon (on land - impossible for glider)
    # Flag 1 (GOOD) = outside island polygon (in the ocean - correct)
    # Flag 9 (MISSING) = coordinates not available
    if 'LATITUDE' in ds and 'LONGITUDE' in ds:
        try:
            # Load La Palma island polygon coordinates from JSON
            polygon_file = 'config/shapefiles/LaPalmaDissolve_coords.json'
            with open(polygon_file, 'r') as f:
                coords_list = json.load(f)
            
            # Create Shapely polygon from coordinates (first polygon in list)
            # Coordinates are in format [[lon, lat], [lon, lat], ...]
            polygon_coords = coords_list[0]
            lapalma_polygon = Polygon(polygon_coords)
            
            # Extract latitude and longitude arrays
            lat_arr = ds['LATITUDE'].values
            lon_arr = ds['LONGITUDE'].values
            
            # Initialize LAND_QC array with flag 9 (missing)
            land_qc = np.full(len(lat_arr), 9, dtype=np.int8)
            
            # Check each point
            for i in range(len(lat_arr)):
                if not (np.isnan(lat_arr[i]) or np.isnan(lon_arr[i])):
                    point = Point(lon_arr[i], lat_arr[i])
                    if lapalma_polygon.contains(point):
                        land_qc[i] = 4  # BAD - inside island (on land)
                    else:
                        land_qc[i] = 1  # GOOD - outside island (in ocean)
            
            # Create dataframe with LAND_QC
            land_qc_df = pd.DataFrame({
                'time': compact['time'],
                'LAND_QC': land_qc
            })
            compact = compact.merge(land_qc_df, on='time', how='left')
            _log.info(f"LAND_QC: Checked {len(lat_arr)} points against La Palma polygon")
            
            # Count flags for logging
            n_good = np.sum(land_qc == 1)
            n_bad = np.sum(land_qc == 4)
            n_missing = np.sum(land_qc == 9)
            _log.info(f"LAND_QC results: {n_good} in ocean (GOOD), {n_bad} on land (BAD), {n_missing} missing")
            
        except FileNotFoundError:
            _log.warning(f"La Palma polygon file not found: {polygon_file}. LAND_QC column will be empty.")
            compact['LAND_QC'] = ''
        except Exception as e:
            _log.warning(f"Error computing LAND_QC: {e}. LAND_QC column will be empty.")
            compact['LAND_QC'] = ''
    else:
        _log.warning("LATITUDE or LONGITUDE not found in dataset. LAND_QC column will be empty.")
        compact['LAND_QC'] = ''

    # merge in requested variable values so QC can be inspected alongside
    value_vars = ['LATITUDE', 'LONGITUDE', 'TEMP', 'CNDC', 'PRES', 'CHLA', 'TURB', 'DOXY', 'PSAL', 'potential_density']
    present = [v for v in value_vars if v in ds]
    
    # Also check for legacy 'temperature' and treat it as 'TEMP'
    if 'temperature' in ds and 'TEMP' not in present:
        present.append('temperature')
    
    # Also check for legacy 'salinity' and treat it as 'PSAL'
    if 'salinity' in ds and 'PSAL' not in present:
        present.append('salinity')
    
    if present:
        # create dataframe of values (may include TIME named differently)
        df_vals = ds[present].to_dataframe().reset_index()
        # ensure time column normalized to 'time'
        vals_time_col = None
        for c in df_vals.columns:
            if c.lower() == 'time':
                vals_time_col = c
                break
        if vals_time_col is None:
            vals_time_col = df_vals.columns[0]
        df_vals = df_vals.rename(columns={vals_time_col: 'time'})
        
        # Rename 'temperature' to 'TEMP' for standardization
        if 'temperature' in df_vals.columns:
            df_vals = df_vals.rename(columns={'temperature': 'TEMP'})
            _log.info("Renamed 'temperature' to 'TEMP' for standardization")
        
        # Rename 'salinity' to 'PSAL' for standardization and convert units
        if 'salinity' in df_vals.columns:
            # Salinity is already in PSU, just rename
            df_vals = df_vals.rename(columns={'salinity': 'PSAL'})
            _log.info("Renamed salinity to PSAL (already in PSU units)")
        
        # Rename 'potential_density' to 'POTDEN' for standardization
        if 'potential_density' in df_vals.columns:
            df_vals = df_vals.rename(columns={'potential_density': 'POTDEN'})
            _log.info("Renamed potential_density to POTDEN for standardization")
        
        # merge values into compact by time
        compact = compact.merge(df_vals, on='time', how='left')

    # ========================================================================
    # SPIKE QC: Detect abrupt changes in T, S, O2
    # ========================================================================
    # Spike test compares each observation V2 with adjacent observations V1 and V3
    # VR = |V2 - (V3 + V1)/2| - |(V3 - V1)/2|
    # Flag 4 (BAD) if VR exceeds threshold (depth-dependent)
    # Flag 0 (NOT EVALUATED) for first and last observations
    # Flag 9 (MISSING) if any of V1, V2, V3, or PRES is NaN
    
    def spike_qc_generic(values, pressure, r1_threshold, r2_threshold, var_name):
        """Generic spike test for any variable.
        
        Args:
            values: numpy array of variable values
            pressure: numpy array of pressure values (dbar)
            r1_threshold: threshold for pressure < 500 dbar
            r2_threshold: threshold for pressure >= 500 dbar
            var_name: variable name for logging
        
        Returns:
            numpy array of QC flags (0, 1, 4, 9)
        """
        n = len(values)
        qc_flags = np.full(n, 0, dtype=int)  # Initialize with 0 (not evaluated)
        
        if n < 3:
            _log.warning(f"Spike QC for {var_name}: dataset has < 3 points, all flagged as 0")
            return qc_flags
        
        # Loop through observations (skip first and last)
        for i in range(1, n - 1):
            V1 = values[i - 1]
            V2 = values[i]
            V3 = values[i + 1]
            P = pressure[i]
            
            # Check for missing data
            if np.isnan(V1) or np.isnan(V2) or np.isnan(V3) or np.isnan(P):
                qc_flags[i] = 0  # NOT EVALUATED (missing data)
                continue
            
            # Calculate VR (spike metric)
            VR = abs(V2 - (V3 + V1) / 2.0) - abs((V3 - V1) / 2.0)
            
            # Determine threshold based on pressure
            if P < 500.0:
                threshold = r1_threshold  # Shallow/intermediate waters
            else:
                threshold = r2_threshold  # Deep waters
            
            # Assign flag
            if VR >= threshold:
                qc_flags[i] = 4  # SPIKE detected
            else:
                qc_flags[i] = 1  # Reasonable
        
        # First and last remain 0 (not evaluated)
        return qc_flags
    
    # TEMP Spike QC
    if 'TEMP' in compact.columns and 'PRES' in compact.columns:
        temp_values = compact['TEMP'].values
        pres_values = compact['PRES'].values
        temp_spike_qc = spike_qc_generic(
            temp_values, pres_values,
            r1_threshold=6.0,   # °C for PRES < 500 dbar
            r2_threshold=2.0,   # °C for PRES >= 500 dbar
            var_name='TEMP'
        )
        compact['TEMP_Spike_QC'] = temp_spike_qc
        _log.info('Applied TEMP Spike QC: R1=6.0°C (P<500), R2=2.0°C (P≥500)')
    else:
        compact['TEMP_Spike_QC'] = 0
        _log.info("TEMP or PRES not available, TEMP_Spike_QC set to 0")
    
    # PSAL Spike QC
    if 'PSAL' in compact.columns and 'PRES' in compact.columns:
        psal_values = compact['PSAL'].values
        pres_values = compact['PRES'].values
        psal_spike_qc = spike_qc_generic(
            psal_values, pres_values,
            r1_threshold=0.9,   # PSU for PRES < 500 dbar
            r2_threshold=0.3,   # PSU for PRES >= 500 dbar
            var_name='PSAL'
        )
        compact['PSAL_Spike_QC'] = psal_spike_qc
        _log.info('Applied PSAL Spike QC: R1=0.9 PSU (P<500), R2=0.3 PSU (P≥500)')
    else:
        compact['PSAL_Spike_QC'] = 0
        _log.info("PSAL or PRES not available, PSAL_Spike_QC set to 0")
    
    # DOXY Spike QC
    if 'DOXY' in compact.columns and 'PRES' in compact.columns:
        doxy_values = compact['DOXY'].values
        pres_values = compact['PRES'].values
        doxy_spike_qc = spike_qc_generic(
            doxy_values, pres_values,
            r1_threshold=50.0,  # µmol/kg for PRES < 500 dbar
            r2_threshold=25.0,  # µmol/kg for PRES >= 500 dbar
            var_name='DOXY'
        )
        compact['DOXY_Spike_QC'] = doxy_spike_qc
        _log.info('Applied DOXY Spike QC: R1=50.0 µmol/kg (P<500), R2=25.0 µmol/kg (P≥500)')
    else:
        compact['DOXY_Spike_QC'] = 0
        _log.info("DOXY or PRES not available, DOXY_Spike_QC set to 0")

    # ============================================================================
    # NEGATIVE SPIKE TEST for CHLA and TURB (5-point median method)
    # ============================================================================
    # This test identifies negative anomalous values using a flexible 5-observation moving window
    # Formula: RES = V2 - median(V0, V1, V2, V3, V4)
    # Flag 4 (SPIKE) if RES < 2 * percentile_10(RES)
    # Flag 0 (NOT EVALUATED) for first and last 2 observations
    # Flag 1 (GOOD) otherwise
    # Flag 9 (MISSING) if V2 itself is NaN or cannot find enough valid neighbors
    # NOTE: If neighbors are NaN, the algorithm searches further back/forward to find valid values
    
    def spike_qc_negative_5point(values, var_name):
        """Negative spike test using flexible 5-point median for CHLA and TURB.
        
        This function uses a flexible window: if a neighbor value is NaN, it searches
        further back or forward to find valid values for the median calculation.
        
        Args:
            values: numpy array of variable values
            var_name: variable name for logging
        
        Returns:
            numpy array of QC flags (0, 1, 4, 9)
        """
        n = len(values)
        qc_flags = np.full(n, 0, dtype=int)  # Initialize with 0 (not evaluated)
        
        if n < 5:
            _log.warning(f"Negative spike QC for {var_name}: dataset has < 5 points, all flagged as 0")
            return qc_flags
        
        def find_valid_neighbors(idx, values, n_before=2, n_after=2):
            """Find valid (non-NaN) neighbors around index idx.
            
            Args:
                idx: current index
                values: array of values
                n_before: number of valid values needed before idx
                n_after: number of valid values needed after idx
            
            Returns:
                tuple: (list of values before, list of values after) or (None, None) if not enough found
            """
            before_vals = []
            after_vals = []
            
            # Search backwards for n_before valid values
            search_idx = idx - 1
            while len(before_vals) < n_before and search_idx >= 0:
                if not np.isnan(values[search_idx]):
                    before_vals.insert(0, values[search_idx])  # Insert at beginning to maintain order
                search_idx -= 1
            
            # Search forwards for n_after valid values
            search_idx = idx + 1
            while len(after_vals) < n_after and search_idx < len(values):
                if not np.isnan(values[search_idx]):
                    after_vals.append(values[search_idx])
                search_idx += 1
            
            # Return None if we couldn't find enough valid values
            if len(before_vals) < n_before or len(after_vals) < n_after:
                return None, None
            
            return before_vals, after_vals
        
        # Step 1: Calculate RES for all valid observations
        res_list = []
        res_indices = []
        
        for i in range(2, n - 2):  # Skip first 2 and last 2
            V2 = values[i]  # Current value
            
            # If V2 itself is missing, mark as not evaluated (0)
            if np.isnan(V2):
                qc_flags[i] = 0  # NOT EVALUATED (missing data)
                continue
            
            # Find valid neighbors (2 before, 2 after)
            before_vals, after_vals = find_valid_neighbors(i, values, n_before=2, n_after=2)
            
            # If we couldn't find enough valid neighbors, mark as not evaluated (0)
            if before_vals is None or after_vals is None:
                qc_flags[i] = 0  # NOT EVALUATED (not enough valid neighbors)
                continue
            
            # Calculate RES = V2 - median(before_vals + [V2] + after_vals)
            window_values = before_vals + [V2] + after_vals
            median_val = np.median(window_values)
            res = V2 - median_val
            res_list.append(res)
            res_indices.append(i)
        
        # Step 2: Calculate threshold (2 * 10th percentile of RES)
        if len(res_list) == 0:
            _log.warning(f"Negative spike QC for {var_name}: no valid RES values calculated")
            return qc_flags
        
        percentile_10 = np.percentile(res_list, 10)
        threshold = 2.0 * percentile_10
        
        _log.info(f"{var_name} Negative Spike QC: P10(RES)={percentile_10:.4f}, threshold={threshold:.4f}")
        
        # Step 3: Assign flags based on threshold
        for idx, res in zip(res_indices, res_list):
            if res < threshold:
                qc_flags[idx] = 4  # SPIKE (negative anomaly)
            else:
                qc_flags[idx] = 1  # GOOD
        
        # First 2 and last 2 remain 0 (not evaluated)
        return qc_flags
    
    # CHLA Negative Spike QC
    if 'CHLA' in compact.columns:
        chla_values = compact['CHLA'].values
        chla_spike_qc = spike_qc_negative_5point(chla_values, var_name='CHLA')
        compact['CHLA_Spike_QC'] = chla_spike_qc
        
        # Log statistics
        n_spikes = np.sum(chla_spike_qc == 4)
        n_good = np.sum(chla_spike_qc == 1)
        n_missing = np.sum(chla_spike_qc == 9)
        _log.info(f'Applied CHLA Negative Spike QC: {n_spikes} spikes, {n_good} good, {n_missing} missing')
    else:
        compact['CHLA_Spike_QC'] = 0
        _log.info("CHLA not available, CHLA_Spike_QC set to 0")
    
    # TURB Negative Spike QC
    if 'TURB' in compact.columns:
        turb_values = compact['TURB'].values
        turb_spike_qc = spike_qc_negative_5point(turb_values, var_name='TURB')
        compact['TURB_Spike_QC'] = turb_spike_qc
        
        # Log statistics
        n_spikes = np.sum(turb_spike_qc == 4)
        n_good = np.sum(turb_spike_qc == 1)
        n_missing = np.sum(turb_spike_qc == 9)
        _log.info(f'Applied TURB Negative Spike QC: {n_spikes} spikes, {n_good} good, {n_missing} missing')
    else:
        compact['TURB_Spike_QC'] = 0
        _log.info("TURB not available, TURB_Spike_QC set to 0")

    # ============================================================================
    # VERTICAL GRADIENT TEST for DOXY (3-point method with flexible neighbors)
    # ============================================================================
    # This test identifies excessive vertical variability in dissolved oxygen
    # Formula: VR = |V2 - (V3+V1)/2|
    # Flag 4 if VR ≥ 50 µmol/kg (P<500 dbar) or VR ≥ 25 µmol/kg (P≥500 dbar)
    # Complementary to spike test: detects rapid vertical transitions (not just isolated spikes)
    
    def doxy_gradient_qc_3point(doxy_values, pressure_values):
        """
        Vertical gradient test for DOXY using flexible 3-point window
        
        Detects excessive vertical changes in dissolved oxygen that may indicate:
        - Sensor calibration drift
        - Unnaturally rapid vertical transitions
        - Sensor response time issues
        
        Args:
            doxy_values: DOXY array (µmol/kg)
            pressure_values: PRES array (dbar)
        
        Returns:
            QC flags: 0=not evaluated, 1=good, 4=excessive gradient, 9=missing
        """
        n = len(doxy_values)
        qc_flags = np.zeros(n, dtype=int)
        
        if n < 3:
            _log.warning("DOXY Gradient QC: dataset has < 3 points, all flagged as 0")
            return qc_flags
        
        def find_valid_neighbor_before(values, idx):
            """Find first valid value before idx"""
            search_idx = idx - 1
            while search_idx >= 0:
                if not np.isnan(values[search_idx]):
                    return values[search_idx]
                search_idx -= 1
            return None
        
        def find_valid_neighbor_after(values, idx):
            """Find first valid value after idx"""
            search_idx = idx + 1
            while search_idx < len(values):
                if not np.isnan(values[search_idx]):
                    return values[search_idx]
                search_idx += 1
            return None
        
        # Loop through all observations except first and last
        for i in range(1, n - 1):
            V2 = doxy_values[i]  # Current value
            P = pressure_values[i]
            
            # If V2 or P is missing, mark as not evaluated (0)
            if np.isnan(V2) or np.isnan(P):
                qc_flags[i] = 0
                continue
            
            # Find valid neighbors (flexible search)
            V1 = find_valid_neighbor_before(doxy_values, i)
            V3 = find_valid_neighbor_after(doxy_values, i)
            
            # If cannot find both neighbors, mark as not evaluated (0)
            if V1 is None or V3 is None:
                qc_flags[i] = 0
                continue
            
            # Calculate VR = |V2 - (V3+V1)/2|
            VR = abs(V2 - (V3 + V1) / 2.0)
            
            # Depth-dependent threshold
            threshold = 50.0 if P < 500.0 else 25.0
            
            # Assign flag
            if VR >= threshold:
                qc_flags[i] = 4  # Excessive vertical gradient
            else:
                qc_flags[i] = 1  # Good
        
        # First and last remain 0 (not evaluated)
        return qc_flags
    
    # Apply DOXY Gradient QC
    if 'DOXY' in compact.columns and 'PRES' in compact.columns:
        doxy_values = compact['DOXY'].values
        pres_values = compact['PRES'].values
        doxy_gradient_qc = doxy_gradient_qc_3point(doxy_values, pres_values)
        compact['DOXY_Gradient_QC'] = doxy_gradient_qc
        
        # Log statistics
        n_grad = np.sum(doxy_gradient_qc == 4)
        n_good = np.sum(doxy_gradient_qc == 1)
        n_missing = np.sum(doxy_gradient_qc == 0)
        _log.info(f'Applied DOXY Gradient QC: {n_grad} excessive gradients, {n_good} good, {n_missing} missing')
        _log.info('DOXY Gradient QC: Thresholds 50 µmol/kg (P<500), 25 µmol/kg (P≥500)')
    else:
        compact['DOXY_Gradient_QC'] = 0
        _log.info("DOXY or PRES not available, DOXY_Gradient_QC set to 0")

    # ------------------------------------------------------------------
    # PRES Max QC
    # Check PRES against 1000 dbar + 10% margin (1100 dbar)
    # PRES missing -> 0, PRES <= 1100 -> 1, PRES > 1100 -> 4
    # ------------------------------------------------------------------
    if 'PRES' in compact.columns:
        pres_vals = compact['PRES'].values
        n = len(pres_vals)
        pres_max_q = np.zeros(n, dtype=int)
        for i in range(n):
            p = pres_vals[i]
            if np.isnan(p):
                pres_max_q[i] = 0
            else:
                pres_max_q[i] = 1 if p <= 1100.0 else 4
        compact['PRES_Max_QC'] = pres_max_q
        _log.info('Applied PRES_Max_QC: PRES<=1100 -> 1, PRES>1100 -> 4')
    else:
        compact['PRES_Max_QC'] = 0
        _log.info('PRES not available, PRES_Max_QC set to 0')

    # ------------------------------------------------------------------
    # DENSITY INVERSION QC (TEOS-10 based)
    # Check for potential density inversions between consecutive observations
    # Uses gsw (Gibbs SeaWater TEOS-10) library to compute potential density
    # Evaluates both downward (bajada) and upward (subida) density gradients
    # Flag 4 if Δ_bajada ≥ 0.03 OR Δ_subida ≤ -0.03; Flag 1 otherwise; Flag 0 if missing/not evaluated
    # ------------------------------------------------------------------
    try:
        import gsw
    except ImportError:
        _log.error("gsw library not found. Please install it: conda install -c conda-forge gsw")
        gsw = None
    
    if gsw and 'TEMP' in compact.columns and 'PSAL' in compact.columns and 'PRES' in compact.columns and 'LATITUDE' in compact.columns and 'LONGITUDE' in compact.columns:
        temp_vals = compact['TEMP'].values
        psal_vals = compact['PSAL'].values
        pres_vals = compact['PRES'].values
        lat_vals = compact['LATITUDE'].values
        lon_vals = compact['LONGITUDE'].values
        n = len(temp_vals)
        
        # Initialize QC flags (0 = not evaluated)
        temp_density_q = np.zeros(n, dtype=int)
        pres_density_q = np.zeros(n, dtype=int)
        psal_density_q = np.zeros(n, dtype=int)
        
        density_inversions = 0
        
        # Loop through consecutive pairs
        for i in range(n - 1):
            temp_i = temp_vals[i]
            temp_i1 = temp_vals[i + 1]
            psal_i = psal_vals[i]
            psal_i1 = psal_vals[i + 1]
            pres_i = pres_vals[i]
            pres_i1 = pres_vals[i + 1]
            lat_i = lat_vals[i]
            lat_i1 = lat_vals[i + 1]
            lon_i = lon_vals[i]
            lon_i1 = lon_vals[i + 1]
            
            # Check for missing data
            if (np.isnan(temp_i) or np.isnan(temp_i1) or 
                np.isnan(psal_i) or np.isnan(psal_i1) or 
                np.isnan(pres_i) or np.isnan(pres_i1) or
                np.isnan(lat_i) or np.isnan(lat_i1) or
                np.isnan(lon_i) or np.isnan(lon_i1)):
                temp_density_q[i] = 0
                pres_density_q[i] = 0
                psal_density_q[i] = 0
                continue
            
            # Calculate reference pressure (mean of two observations)
            pref = (pres_i + pres_i1) / 2.0
            
            try:
                # Calculate Absolute Salinity from Practical Salinity
                sa_i = gsw.SA_from_SP(psal_i, pres_i, lon_i, lat_i)
                sa_i1 = gsw.SA_from_SP(psal_i1, pres_i1, lon_i1, lat_i1)
                
                # Calculate potential density at reference pressure
                rho_i = gsw.pot_rho_t_exact(sa_i, temp_i, pres_i, pref)
                rho_i1 = gsw.pot_rho_t_exact(sa_i1, temp_i1, pres_i1, pref)
                
                # Calculate density differences
                delta_bajada = rho_i - rho_i1  # Downward gradient
                delta_subida = rho_i1 - rho_i  # Upward gradient (= -delta_bajada)
                
                # Check for inversions (threshold = 0.03 kg/m³)
                threshold = 0.03
                if delta_bajada >= threshold or delta_subida <= -threshold:
                    # Density inversion detected
                    temp_density_q[i] = 4
                    pres_density_q[i] = 4
                    psal_density_q[i] = 4
                    density_inversions += 1
                else:
                    # No inversion
                    temp_density_q[i] = 1
                    pres_density_q[i] = 1
                    psal_density_q[i] = 1
            except Exception as e:
                _log.warning(f"Error computing density at index {i}: {e}. Setting QC to 0.")
                temp_density_q[i] = 0
                pres_density_q[i] = 0
                psal_density_q[i] = 0
        
        # Last observation is not evaluated
        temp_density_q[n - 1] = 0
        pres_density_q[n - 1] = 0
        psal_density_q[n - 1] = 0
        
        # Assign to compact dataframe
        compact['TEMP_Density_QC'] = temp_density_q
        compact['PRES_Density_QC'] = pres_density_q
        compact['PSAL_Density_QC'] = psal_density_q
        
        _log.info(f'Applied Density Inversion QC: {density_inversions} inversions detected (threshold: 0.03 kg/m³)')
    else:
        # gsw not available or required variables missing
        if not gsw:
            _log.error("gsw library is required for Density Inversion QC. Install: conda install -c conda-forge gsw")
            compact['TEMP_Density_QC'] = 0
            compact['PRES_Density_QC'] = 0
            compact['PSAL_Density_QC'] = 0
        else:
            missing_vars = []
            if 'TEMP' not in compact.columns:
                missing_vars.append('TEMP')
            if 'PSAL' not in compact.columns:
                missing_vars.append('PSAL')
            if 'PRES' not in compact.columns:
                missing_vars.append('PRES')
            if 'LATITUDE' not in compact.columns:
                missing_vars.append('LATITUDE')
            if 'LONGITUDE' not in compact.columns:
                missing_vars.append('LONGITUDE')
            _log.warning(f"Density Inversion QC: missing variables {missing_vars}. QC columns set to 0.")
            compact['TEMP_Density_QC'] = 0
            compact['PRES_Density_QC'] = 0
            compact['PSAL_Density_QC'] = 0

    # ------------------------------------------------------------------
    # PRES Increasing QC
    # Compare current PRES with previous observation: if PRES changes (increase or decrease) -> 1, else -> 4
    # First observation or missing -> 0
    # ------------------------------------------------------------------
    if 'PRES' in compact.columns:
        pres_vals = compact['PRES'].values
        n = len(pres_vals)
        pres_inc_q = np.zeros(n, dtype=int)
        for i in range(n):
            if i == 0:
                pres_inc_q[i] = 0
                continue
            v1 = pres_vals[i-1]
            v2 = pres_vals[i]
            if np.isnan(v1) or np.isnan(v2):
                pres_inc_q[i] = 0
            else:
                # If pressure changes (either increases or decreases) it's acceptable -> 1
                if v2 > v1 or v2 < v1:
                    pres_inc_q[i] = 1
                else:
                    pres_inc_q[i] = 4
        compact['PRES_Increasing_QC'] = pres_inc_q
        _log.info('Applied PRES_Increasing_QC: change->1, no-change->4, first/missing->0')
    else:
        compact['PRES_Increasing_QC'] = 0
        _log.info('PRES not available, PRES_Increasing_QC set to 0')

    # ------------------------------------------------------------------
    # STUCK VALUES TEST (Perfiles con valores idénticos)
    # Detects if a variable is "stuck" (all identical values) within a profile segment
    # Profile segments are determined by PRES direction: downcast (increasing) or upcast (decreasing)
    # Flag 4 = all non-NaN values identical (sensor stuck), Flag 1 = values vary, Flag 0 = not in valid profile
    # Applied to: TEMP, CNDC, DOXY, CHLA, TURB
    # ------------------------------------------------------------------
    if 'PRES' in compact.columns:
        pres_vals = compact['PRES'].values
        n = len(pres_vals)
        
        # Identify profile segments based on PRES direction
        # segment_id[i] = which profile segment observation i belongs to
        # profile_direction[i] = 'down' (increasing PRES), 'up' (decreasing PRES), or 'none' (isolated)
        segment_id = np.full(n, -1, dtype=int)
        profile_direction = np.full(n, 'none', dtype=object)
        
        current_segment = 0
        current_direction = None
        
        for i in range(n):
            if i == 0:
                segment_id[i] = current_segment
                continue
            
            p_prev = pres_vals[i - 1]
            p_curr = pres_vals[i]
            
            # Skip if pressure is NaN
            if np.isnan(p_prev) or np.isnan(p_curr):
                segment_id[i] = -1
                continue
            
            # Determine direction
            if p_curr > p_prev:
                new_direction = 'down'
            elif p_curr < p_prev:
                new_direction = 'up'
            else:
                # PRES unchanged; stay in same segment
                new_direction = current_direction
            
            # Check if direction changed
            if current_direction is not None and new_direction != current_direction and new_direction != current_direction:
                current_segment += 1
            
            current_direction = new_direction
            segment_id[i] = current_segment
            profile_direction[i] = new_direction
        
        # Now apply stuck values test to each variable
        stuck_vars = ['TEMP', 'CNDC', 'DOXY', 'CHLA', 'TURB']
        stuck_qc_dict = {}
        
        for var in stuck_vars:
            if var not in compact.columns:
                stuck_qc_dict[var] = np.zeros(n, dtype=int)
                continue
            
            var_vals = compact[var].values
            stuck_qc = np.zeros(n, dtype=int)  # Initialize all as 0 (not evaluated)
            
            # Process each segment
            for seg_id in np.unique(segment_id):
                if seg_id < 0:
                    # Invalid segment
                    continue
                
                # Get indices in this segment
                seg_mask = segment_id == seg_id
                seg_indices = np.where(seg_mask)[0]
                
                if len(seg_indices) == 0:
                    continue
                
                # Get values in this segment (non-NaN only)
                seg_values = var_vals[seg_indices]
                valid_mask = ~np.isnan(seg_values)
                valid_values = seg_values[valid_mask]
                
                if len(valid_values) == 0:
                    # All NaN in this segment
                    stuck_qc[seg_indices] = 0
                    continue
                
                # Check if all valid values are identical
                if np.all(valid_values == valid_values[0]):
                    # All values are the same -> stuck
                    stuck_qc[seg_indices] = 4
                else:
                    # Values vary -> good
                    stuck_qc[seg_indices] = 1
            
            stuck_qc_dict[var] = stuck_qc
        
        # Assign to compact dataframe
        compact['TEMP_Stuck_QC'] = stuck_qc_dict['TEMP']
        compact['CNDC_Stuck_QC'] = stuck_qc_dict['CNDC']
        compact['DOXY_Stuck_QC'] = stuck_qc_dict['DOXY']
        compact['CHLA_Stuck_QC'] = stuck_qc_dict['CHLA']
        compact['TURB_Stuck_QC'] = stuck_qc_dict['TURB']
        
        # Log statistics
        for var in stuck_vars:
            stuck_count = np.sum(stuck_qc_dict[var] == 4)
            if stuck_count > 0:
                _log.info(f'Stuck Values Test for {var}: {stuck_count} observations with identical values detected')
        
        _log.info('Applied Stuck Values Test based on profile segments (down/upcast)')
    else:
        compact['TEMP_Stuck_QC'] = 0
        compact['CNDC_Stuck_QC'] = 0
        compact['DOXY_Stuck_QC'] = 0
        compact['CHLA_Stuck_QC'] = 0
        compact['TURB_Stuck_QC'] = 0
        _log.info('PRES not available, Stuck Values QC set to 0')

    # ------------------------------------------------------------------
    # SURFACE QC for CHLA and TURB
    # If PRES <= 5 dbar -> flag 4 (bad at surface), else flag 1 (reasonable)
    # If PRES missing or variable missing -> flag 0 (not evaluated)
    # ------------------------------------------------------------------
    # CHLA Surface QC
    if 'CHLA' in compact.columns and 'PRES' in compact.columns:
        pres_vals = compact['PRES'].values
        chla_vals = compact['CHLA'].values
        n = len(pres_vals)
        chla_surface_q = np.zeros(n, dtype=int)
        for i in range(n):
            p = pres_vals[i]
            # not evaluated if pressure or value missing
            if np.isnan(p) or (np.isnan(chla_vals[i]) if not np.isscalar(chla_vals[i]) else False):
                chla_surface_q[i] = 0
            else:
                chla_surface_q[i] = 4 if p <= 5.0 else 1
        compact['CHLA_Surface_QC'] = chla_surface_q
        _log.info('Applied CHLA Surface QC: PRES<=5 -> 4, PRES>5 -> 1')
    else:
        compact['CHLA_Surface_QC'] = 0
        _log.info('CHLA or PRES not available, CHLA_Surface_QC set to 0')

    # TURB Surface QC
    if 'TURB' in compact.columns and 'PRES' in compact.columns:
        pres_vals = compact['PRES'].values
        turb_vals = compact['TURB'].values
        n = len(pres_vals)
        turb_surface_q = np.zeros(n, dtype=int)
        for i in range(n):
            p = pres_vals[i]
            if np.isnan(p) or (np.isnan(turb_vals[i]) if not np.isscalar(turb_vals[i]) else False):
                turb_surface_q[i] = 0
            else:
                turb_surface_q[i] = 4 if p <= 5.0 else 1
        compact['TURB_Surface_QC'] = turb_surface_q
        _log.info('Applied TURB Surface QC: PRES<=5 -> 4, PRES>5 -> 1')
    else:
        compact['TURB_Surface_QC'] = 0
        _log.info('TURB or PRES not available, TURB_Surface_QC set to 0')

    # For column ordering, replace 'temperature' and 'salinity' with standardized names
    present_standardized = []
    for v in present:
        if v == 'salinity':
            present_standardized.append('PSAL')
        elif v == 'temperature':
            present_standardized.append('TEMP')
        elif v == 'potential_density':
            present_standardized.append('POTDEN')
        else:
            present_standardized.append(v)
    
    # Create Na_QC (missing value) columns for specified variables
    # Flag 1 = value present, Flag 9 = value missing (NaN)
    na_qc_vars = ['TEMP', 'CNDC', 'DOXY', 'CHLA', 'TURB']
    for var in na_qc_vars:
        na_qc_col = f'{var}_Na_QC'
        if var in compact.columns:
            # Check if value is present (1) or missing (9)
            compact[na_qc_col] = compact[var].notna().astype(int).replace({1: 1, 0: 9})
        else:
            # If variable doesn't exist, leave Na_QC column empty
            compact[na_qc_col] = ''
    
    # reorder columns to TIME, DEPTH, requested values, then Range_QC columns, then Na_QC columns, then Sensor_QC, then LAND_QC, then Spike_QC, then Gradient_QC
    range_qc_cols = ['Date_QC', 'Location_QC', 'TEMP_Range_QC', 'CNDC_Range_QC', 'CHLA_Range_QC', 'TURB_Range_QC', 'DOXY_Range_QC', 'PSAL_Range_QC']
    na_qc_cols = ['TEMP_Na_QC', 'CNDC_Na_QC', 'DOXY_Na_QC', 'CHLA_Na_QC', 'TURB_Na_QC']
    sensor_qc_cols = ['TEMP_Sensor_QC', 'CNDC_Sensor_QC', 'DOXY_Sensor_QC', 'CHLA_Sensor_QC', 'TURB_Sensor_QC']
    spike_qc_cols = ['TEMP_Spike_QC', 'PSAL_Spike_QC', 'DOXY_Spike_QC', 'CHLA_Spike_QC', 'TURB_Spike_QC']
    land_qc_cols = ['LAND_QC']
    gradient_qc_cols = ['DOXY_Gradient_QC']
    surface_qc_cols = ['CHLA_Surface_QC', 'TURB_Surface_QC']
    pres_qc_cols = ['PRES_Max_QC']
    density_qc_cols = ['TEMP_Density_QC', 'PRES_Density_QC', 'PSAL_Density_QC']
    increasing_qc_cols = ['PRES_Increasing_QC']
    stuck_qc_cols = ['TEMP_Stuck_QC', 'CNDC_Stuck_QC', 'DOXY_Stuck_QC', 'CHLA_Stuck_QC', 'TURB_Stuck_QC']
    # Order: ... spike QC, surface QC, gradient QC, PRES_Max_QC, Density QC, PRES_Increasing_QC, then Stuck Values QC
    cols = ['time', 'depth'] + present_standardized + range_qc_cols + na_qc_cols + sensor_qc_cols + land_qc_cols + spike_qc_cols + surface_qc_cols + gradient_qc_cols + pres_qc_cols + density_qc_cols + increasing_qc_cols + stuck_qc_cols
    # keep only columns that exist in compact
    cols = [c for c in cols if c in compact.columns]
    compact = compact[cols]

    # Convert all QC flag columns: replace NaN (and empty strings) with 0 for QC-only columns
    # We target only the QC columns (range, sensor, spike, gradient, date/location, LAND_QC)
    qc_columns = range_qc_cols + sensor_qc_cols + spike_qc_cols + surface_qc_cols + gradient_qc_cols + pres_qc_cols + density_qc_cols + increasing_qc_cols + stuck_qc_cols + land_qc_cols + ['Date_QC', 'Location_QC']
    for col in qc_columns:
        if col in compact.columns:
            # Normalize empty strings to NaN, then fill NaN with 0 (not evaluated/missing)
            compact[col] = compact[col].replace('', np.nan)
            compact[col] = pd.to_numeric(compact[col], errors='coerce').fillna(0).astype(int)

    # Uppercase TIME and DEPTH column headers as requested
    if 'time' in compact.columns:
        compact = compact.rename(columns={'time': 'TIME'})
    if 'depth' in compact.columns:
        compact = compact.rename(columns={'depth': 'DEPTH'})

    if csv_path is None:
        # Use a canonical variables-level filename rather than per-variable by default
        csv_path = f'output/analysis/seaexplorer_qc_variables.csv'
    import os
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    compact.to_csv(csv_path, index=False, na_rep='NaN')
    _log.info('Wrote compact QC CSV to %s', csv_path)
    return csv_path


def export_full_qc_csv(ds: xr.Dataset, varname: Optional[str] = None, csv_path: str = None):
    """Export a fuller CSV containing time, the variable value, QC flag, and coords if present.

    Columns produced (when available): time, <varname>, QC/<Varname>, latitude, longitude, depth
    """
    # Require standardized 'TEMP' for full CSV
    if varname is None:
        if 'TEMP' in ds:
            varname = 'TEMP'
            _log.info("Using standardized temperature variable '%s' for full CSV", varname)
        else:
            raise KeyError("Variable 'TEMP' not found in dataset; full CSV export requires standardized 'TEMP'")
    qc_name = f"{varname}_qc"
    if varname not in ds:
        raise KeyError(f"Variable '{varname}' not found in dataset")
    if qc_name not in ds:
        raise KeyError(f"QC variable '{qc_name}' not found in dataset. Run range_qc_temperature first")

    # Build dataframe from variable and qc, preserving time index
    df_var = ds[varname].to_dataframe(name=varname).reset_index()
    display_name = varname if varname.isupper() else varname.capitalize()
    df_qc = ds[qc_name].to_dataframe(name=f'QC/{display_name}').reset_index()
    # Merge on time (and any other indices)
    # prefer simple merge on 'time' when present
    if 'time' in df_var.columns and 'time' in df_qc.columns:
        df = df_var.merge(df_qc, on='time')
    else:
        # fallback: concatenate side-by-side
        df = df_var.join(df_qc.set_index(df_qc.columns[0]), how='left')

    # add coords if present in the dataset
    for coord in ('latitude', 'longitude', 'depth'):
        if coord in ds.coords and coord not in df.columns:
            # write coord values aligned by index
            df[coord] = ds[coord].values

    if csv_path is None:
        csv_path = f'output/analysis/seaexplorer_qc_variables_full.csv'
    import os
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    df.to_csv(csv_path, index=False)
    _log.info('Wrote full QC CSV to %s', csv_path)
    return csv_path


def _cli_apply_to_netcdf(nc_in: str, nc_out: str = None):
    """Simple CLI helper: open nc_in, apply temp QC.

    By default this writes only a separate CSV with QC flags and does NOT
    overwrite the input NetCDF. If `nc_out` is provided, the NetCDF will be
    overwritten with the QC variable added.
    """
    ds = xr.open_dataset(nc_in)
    # Require TEMP-only behavior
    if 'TEMP' not in ds:
        raise SystemExit("Input NetCDF does not contain standardized 'TEMP' variable required by this QC")
    ds = range_qc_temperature(ds, varname='TEMP')
    # Default CSV path (project structure) - use the canonical variables filename
    default_csv = f'output/analysis/seaexplorer_qc_variables.csv'
    # Ensure output directory exists
    import os
    os.makedirs(os.path.dirname(default_csv), exist_ok=True)
    # write CSV by default
    csv_path = nc_out if (nc_out and nc_out.lower().endswith('.csv')) else default_csv
    # No extra vars supplied by the simple CLI helper; callers may use the
    # module-level CLI to pass --extra if desired.
    export_qc_to_csv(ds, varname='TEMP', csv_path=csv_path)
    print(f'Wrote QC CSV to {csv_path}')
    # Overwrite NetCDF only if nc_out is provided and endswith .nc
    if nc_out and nc_out.lower().endswith('.nc'):
        tmp = nc_out + '.tmp'
        ds.to_netcdf(tmp)
        os.replace(tmp, nc_out)
        print(f'Overwrote NetCDF with QC variable: {nc_out}')


if __name__ == '__main__':
    import argparse

    p = argparse.ArgumentParser(description='Apply temperature range QC to NetCDF file(s)')
    p.add_argument('ncfile', nargs='?', help='Input NetCDF file (optional if --latest or --all used)')
    p.add_argument('--out', help='Optional: write modified NetCDF (path ending with .nc) or CSV (path ending with .csv). If omitted a CSV will be written to output/analysis/seaexplorer_qc_variables.csv')
    p.add_argument('--latest', action='store_true', help='Auto-discover the latest L0 NetCDF in output/l0_data/timeseries/ and process it')
    p.add_argument('--all', action='store_true', help='Process all NetCDF files in output/l0_data/timeseries/ and create per-file CSVs')
    args = p.parse_args()

    # If --latest or --all requested, ignore positional ncfile
    if args.latest or args.all:
        import glob, os
        # Prefer standardized NetCDFs produced by MASTER in output/analysis/
        analysis_pattern = 'output/analysis/*_standard_names.nc'
        analysis_files = sorted(glob.glob(analysis_pattern))
        if analysis_files:
            files = analysis_files
        else:
            # fallback to L0 timeseries folder
            folder = 'output/l0_data/timeseries/'
            pattern = folder + '*.nc'
            files = sorted(glob.glob(pattern))

        if not files:
            print('No NetCDF files found in output/analysis/ or output/l0_data/timeseries/')
            raise SystemExit(1)

        if args.latest:
            target = files[-1]
            _cli_apply_to_netcdf(target, args.out)
        else:
            # process all files; create per-file CSVs in output/analysis/
            for f in files:
                base = os.path.splitext(os.path.basename(f))[0]
                csv_out = f'output/analysis/{base}_qc_variables.csv'
                _cli_apply_to_netcdf(f, csv_out)
    else:

        if not args.ncfile:
            p.print_usage()
            print('\nError: must provide a NetCDF file or use --latest/--all')
            raise SystemExit(2)
        _cli_apply_to_netcdf(args.ncfile, args.out)
