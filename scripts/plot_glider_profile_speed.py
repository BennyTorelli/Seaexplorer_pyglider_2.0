"""Plot glider profile/trajectory colored by speed.

Reads a CSV or NetCDF with time, latitude, longitude and depth/pressure and
creates a figure with:
 - left: depth vs time (depth positive down, inverted y-axis), colored by horizontal speed (m/s)
 - right: depth vs along-track distance (m), colored by horizontal speed (m/s)

Usage:
    python scripts/plot_glider_profile_speed.py --input test_without_nan.csv --output output/analysis/glider_profile_speed.png

If input is NetCDF, the script will try to find variables named 'TIME'/'time',
'LATITUDE'/'latitude', 'LONGITUDE'/'longitude', and 'DEPTH' or 'pressure'.
If only 'pressure' is present, depth(m) is approximated as pressure(dbar).

"""
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import radians, cos, sin, asin, sqrt


def haversine(lon1, lat1, lon2, lat2):
    """Return distance in meters between two lon/lat points."""
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))
    R = 6371000.0  # Earth radius in meters
    return R * c


def compute_speed_from_positions(lat, lon, times):
    """Compute horizontal speed (m/s) between consecutive positions.
    Returns array of speed with same length (first element NaN).
    """
    n = len(lat)
    speed = np.full(n, np.nan)
    if n < 2:
        return speed
    # compute distances and time deltas
    dists = np.zeros(n)
    dt = np.zeros(n)
    for i in range(1, n):
        dists[i] = haversine(lon[i-1], lat[i-1], lon[i], lat[i])
        dt[i] = (times[i] - times[i-1]).total_seconds()
        if dt[i] > 0:
            speed[i] = dists[i] / dt[i]
        else:
            speed[i] = np.nan
    return speed, dists


def load_data(input_path):
    """Load data from CSV or NetCDF into a pandas DataFrame with columns:
    'time' (datetime), 'latitude', 'longitude', 'depth' (meters)
    """
    lower = input_path.lower()
    if lower.endswith('.csv'):
        df = pd.read_csv(input_path, parse_dates=['time'])
        # Normalize column names
        cols = {c.lower(): c for c in df.columns}
        # rename to known names
        mapping = {}
        if 'time' in cols:
            mapping[cols['time']] = 'time'
        if 'latitude' in cols:
            mapping[cols['latitude']] = 'latitude'
        if 'longitude' in cols:
            mapping[cols['longitude']] = 'longitude'
        if 'depth' in cols:
            mapping[cols['depth']] = 'depth'
        if 'pressure' in cols and 'depth' not in mapping:
            mapping[cols['pressure']] = 'pressure'
        df = df.rename(columns=mapping)
        # if pressure present but not depth, approximate depth = pressure
        if 'depth' not in df.columns and 'pressure' in df.columns:
            df['depth'] = df['pressure']
        # ensure required columns
        for req in ('time', 'latitude', 'longitude', 'depth'):
            if req not in df.columns:
                raise ValueError(f"Required column '{req}' not found in {input_path}")
        df = df.sort_values('time').reset_index(drop=True)
        return df[['time', 'latitude', 'longitude', 'depth']]

    else:
        # try xarray for netcdf
        import xarray as xr
        ds = xr.open_dataset(input_path)
        # find names ignoring case
        def find_var(ds, candidates):
            for c in candidates:
                if c in ds:
                    return c
                # case-insensitive
                for v in ds.variables:
                    if v.lower() == c.lower():
                        return v
            return None

        time_var = find_var(ds, ['TIME', 'time'])
        lat_var = find_var(ds, ['LATITUDE', 'latitude', 'lat'])
        lon_var = find_var(ds, ['LONGITUDE', 'longitude', 'lon'])
        depth_var = find_var(ds, ['DEPTH', 'depth', 'pressure', 'PRES', 'pressure'])
        if not all([time_var, lat_var, lon_var, depth_var]):
            raise ValueError('Could not find required vars in NetCDF. Need time, latitude, longitude and depth/pressure')
        times = pd.to_datetime(ds[time_var].values)
        lat = ds[lat_var].values
        lon = ds[lon_var].values
        depth = ds[depth_var].values
        df = pd.DataFrame({'time': times, 'latitude': lat, 'longitude': lon, 'depth': depth})
        df = df.sort_values('time').reset_index(drop=True)
        return df


def plot_profile(df, output_path=None, show=False):
    """Create and save figure."""
    # compute speed
    times = list(df['time'])
    lat = df['latitude'].values.astype(float)
    lon = df['longitude'].values.astype(float)
    depth = df['depth'].values.astype(float)

    speed, dists = compute_speed_from_positions(lat, lon, times)
    cumdist = np.cumsum(np.nan_to_num(dists))

    # Prepare plot
    fig, axs = plt.subplots(1, 2, figsize=(14,6), sharey=True)

    sc = axs[0].scatter(df['time'], depth, c=speed, cmap='viridis', s=20, edgecolors='none')
    axs[0].invert_yaxis()
    axs[0].set_xlabel('Time')
    axs[0].set_ylabel('Depth (m)')
    axs[0].set_title('Depth vs Time colored by horizontal speed (m/s)')
    cb = fig.colorbar(sc, ax=axs[0], label='speed (m/s)')

    sc2 = axs[1].scatter(cumdist, depth, c=speed, cmap='viridis', s=20, edgecolors='none')
    axs[1].invert_yaxis()
    axs[1].set_xlabel('Along-track distance (m)')
    axs[1].set_title('Depth vs Along-track distance colored by speed')
    fig.colorbar(sc2, ax=axs[1], label='speed (m/s)')

    fig.tight_layout()

    if output_path:
        outdir = os.path.dirname(output_path)
        if outdir and not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
        fig.savefig(output_path, dpi=200)
        print(f"Saved figure to {output_path}")
    if show:
        plt.show()
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot glider depth profile colored by speed')
    parser.add_argument('--input', '-i', required=True, help='Input CSV or NetCDF path')
    parser.add_argument('--output', '-o', default='output/analysis/glider_profile_speed.png', help='Output image path')
    parser.add_argument('--show', action='store_true', help='Show the figure interactively')
    args = parser.parse_args()

    df = load_data(args.input)
    plot_profile(df, args.output, show=args.show)
