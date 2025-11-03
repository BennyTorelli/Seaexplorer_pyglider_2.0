"""Backward-compatible shim for scripts.qc_temperature

This module re-exports the implementation from `scripts.qc_variables` so existing
imports keep working while we finish the rename. Once you're happy with the new
module name we can remove this shim.
"""
from .qc_variables import *  # noqa: F401,F403

__all__ = [
    'range_qc_temperature',
    'range_qc_variable',
    'export_qc_to_csv',
    'export_full_qc_csv',
]
"""Simple QC utilities for temperature.

This module prefers the standardized variable name `TEMP` when present in
NetCDF files produced by the MASTER pipeline, but will fall back to the
legacy `temperature` name when necessary. The CLI and helper functions will
auto-detect the available name if the caller doesn't provide one.
"""
from typing import Tuple, Optional
import xarray as xr
import numpy as np
import logging

_log = logging.getLogger(__name__)


def range_qc_temperature(ds: xr.Dataset, varname: Optional[str] = None,
                         valid_min: float = -2.5, valid_max: float = 40.0) -> xr.Dataset:
    """Apply range QC to a temperature variable in an xarray.Dataset.

        Adds or updates a companion quality flag variable named `<varname>_qc` with
        integer flags following the requested convention:
            1 = PASS (in-range, "ottimo")
            4 = FAIL (out-of-range, "non va bene")
            9 = MISSING (nan)

    The function preserves attributes on the original variable and adds
    reasonable attributes to the QC variable.

    Args:
        ds: xarray.Dataset containing the variable to check.
        varname: variable name to QC (default: 'temperature').
        valid_min: minimum valid value (inclusive).
        valid_max: maximum valid value (inclusive).

    Returns:
        The dataset with an added `<varname>_qc` variable (int8).
    """
    # The QC must operate on the standardized variable 'TEMP' only.
    # If varname is provided, honor it only if it equals 'TEMP'.
    if varname is None:
        if 'TEMP' in ds:
            varname = 'TEMP'
            _log.info("Using standardized temperature variable '%s'", varname)
        else:
            raise KeyError("Variable 'TEMP' not found in dataset; this QC operates only on the standardized 'TEMP' variable")
    else:
        if varname != 'TEMP':
            raise KeyError("This QC function only accepts the standardized variable name 'TEMP'")

    data = ds[varname].values
    # initialize flags as MISSING (9)
    flags = np.full(data.shape, 9, dtype=np.int8)

    # Where values are finite, mark PASS or FAIL
    finite_mask = np.isfinite(data)
    if finite_mask.any():
        vals = data[finite_mask]
        pass_mask = (vals >= valid_min) & (vals <= valid_max)
        # PASS -> 1, FAIL -> 4
        flags[finite_mask] = np.where(pass_mask, 1, 4).astype(np.int8)

    qc_name = f"{varname}_qc"
    attrs = {
        'long_name': f'{varname} quality flag',
        'flag_values': np.array([1, 4, 9], dtype=np.int8),
        'flag_meanings': 'pass fail missing',
        'units': '1',
        'method': 'range check',
        'valid_min': 1,
        'valid_max': 9,
    }

    # Insert or replace qc variable
    ds[qc_name] = (ds[varname].dims, flags, attrs)
    _log.info('Applied range QC to %s (%g..%g)', varname, valid_min, valid_max)
    return ds


def range_qc_variable(ds: xr.Dataset, varname: str,
                      valid_min: float, valid_max: float) -> xr.Dataset:
    """Generic range QC for any existing variable name (expects standardized name).

    Adds `<varname>_qc` to the dataset with flags 1=pass,4=fail,9=missing.
    """
    if varname not in ds:
        raise KeyError(f"Variable '{varname}' not found in dataset")

    data = ds[varname].values
    flags = np.full(data.shape, 9, dtype=np.int8)
    finite_mask = np.isfinite(data)
    if finite_mask.any():
        vals = data[finite_mask]
        pass_mask = (vals >= valid_min) & (vals <= valid_max)
        flags[finite_mask] = np.where(pass_mask, 1, 4).astype(np.int8)

    qc_name = f"{varname}_qc"
    attrs = {
        'long_name': f'{varname} quality flag',
        'flag_values': np.array([1, 4, 9], dtype=np.int8),
        'flag_meanings': 'pass fail missing',
        'units': '1',
        'method': 'range check',
        'valid_min': 1,
        'valid_max': 9,
    }
    ds[qc_name] = (ds[varname].dims, flags, attrs)
    _log.info('Applied range QC to %s (%g..%g)', varname, valid_min, valid_max)
    return ds


def export_qc_to_csv(ds: xr.Dataset, varname: str = None, csv_path: str = None):
    """Export a compact QC CSV containing only time, depth, and QC/temperature.

    The output columns will be exactly: time, depth, QC/temperature
    where QC/temperature uses the lowercase legacy label you requested.
    The function auto-detects `TEMP` vs `temperature` when `varname` is None.
    """
    # Require the standardized 'TEMP' variable
    if varname is None:
        if 'TEMP' in ds:
            varname = 'TEMP'
            _log.info("Using standardized temperature variable '%s' for CSV export", varname)
        else:
            raise KeyError("Variable 'TEMP' not found in dataset; export requires standardized 'TEMP'")
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
    compact = compact.rename(columns={time_col: 'time', 'qc': 'QC/TEMP'})
    # attach depth if available (align by index length)
    if depth_vals is not None:
        # depth might be a per-sample array aligned with the dataset
        compact['depth'] = depth_vals
    else:
        compact['depth'] = ''

    # Also compute CNDC QC if CNDC present
    if 'CNDC' in ds:
        ds = range_qc_variable(ds, 'CNDC', 0.0, 7.0)
        cndc_q = ds['CNDC_qc'].to_dataframe(name='cndc_q').reset_index()
        # align on time column
        if time_col in cndc_q.columns:
            cndc_q = cndc_q[[time_col, 'cndc_q']]
            cndc_q = cndc_q.rename(columns={time_col: 'time', 'cndc_q': 'QC/CNDC'})
            compact = compact.merge(cndc_q, on='time', how='left')
        else:
            compact['QC/CNDC'] = ''
    else:
        compact['QC/CNDC'] = ''

    # merge in requested variable values so QC can be inspected alongside
    value_vars = ['LATITUDE', 'LONGITUDE', 'TEMP', 'CNDC', 'PRES', 'CHLA', 'TURB', 'DOXY']
    present = [v for v in value_vars if v in ds]
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
        # merge values into compact by time
        compact = compact.merge(df_vals, on='time', how='left')

    # reorder columns to TIME, DEPTH, requested values, then QC columns
    cols = ['time', 'depth'] + present + ['QC/TEMP', 'QC/CNDC']
    # keep only columns that exist in compact
    cols = [c for c in cols if c in compact.columns]
    compact = compact[cols]

    # Uppercase TIME and DEPTH column headers as requested
    if 'time' in compact.columns:
        compact = compact.rename(columns={'time': 'TIME'})
    if 'depth' in compact.columns:
        compact = compact.rename(columns={'depth': 'DEPTH'})

    if csv_path is None:
        csv_path = f'output/analysis/seaexplorer_qc_{varname}.csv'
    import os
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    compact.to_csv(csv_path, index=False)
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
        csv_path = f'output/analysis/seaexplorer_qc_{varname}_full.csv'
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
    # Default CSV path (project structure)
    default_csv = f'output/analysis/seaexplorer_qc_TEMP.csv'
    # Ensure output directory exists
    import os
    os.makedirs(os.path.dirname(default_csv), exist_ok=True)
    # write CSV by default
    csv_path = nc_out if (nc_out and nc_out.lower().endswith('.csv')) else default_csv
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
    p.add_argument('--out', help='Optional: write modified NetCDF (path ending with .nc) or CSV (path ending with .csv). If omitted a CSV will be written to output/analysis/seaexplorer_qc_temperature.csv')
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
                csv_out = f'output/analysis/{base}_qc_temperature.csv'
                _cli_apply_to_netcdf(f, csv_out)
    else:
        if not args.ncfile:
            p.print_usage()
            print('\nError: must provide a NetCDF file or use --latest/--all')
            raise SystemExit(2)
        _cli_apply_to_netcdf(args.ncfile, args.out)
