"""Standalone runner to apply QC to the latest MASTER-produced NetCDF.

This script does not modify MASTER_pyglider_pipeline.py. It finds the
latest standardized NetCDF in output/analysis/ (pattern '*_standard_names.nc')
or falls back to output/l0_data/timeseries/, applies QC from
`scripts.qc_variables`, and writes the canonical CSV to
`output/analysis/seaexplorer_qc_variables.csv`.
"""
import glob
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.qc_variables import range_qc_temperature, export_qc_to_csv


def find_latest_l0():
    pattern = 'output/analysis/*_standard_names.nc'
    files = sorted(glob.glob(pattern))
    if files:
        return files[-1]
    pattern2 = 'output/l0_data/timeseries/*.nc'
    files2 = sorted(glob.glob(pattern2))
    return files2[-1] if files2 else None


def main():
    nc = find_latest_l0()
    if not nc:
        print('No L0 NetCDF found in output/analysis/ or output/l0_data/timeseries/')
        raise SystemExit(1)

    print(f'Applying QC to: {nc}')
    ds = None
    import xarray as xr
    try:
        ds = xr.open_dataset(nc)
        ds = range_qc_temperature(ds)
        out = export_qc_to_csv(ds, varname='TEMP')
        print(f'Wrote QC CSV: {out}')
    finally:
        if ds is not None:
            ds.close()


if __name__ == '__main__':
    main()
