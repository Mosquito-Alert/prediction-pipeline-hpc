import argparse
from datetime import datetime
import numpy as np
import xarray as xr

# TODO: check why there's more than 1 day.

def main(input_file: str, output_file: str):
    ds = xr.open_dataset(input_file)

    # Calculate mean of d2m and t2m
    d2m_mean = ds['d2m'].mean(dim=['valid_time']) - 273.15   # Convert to Celsius
    t2m_mean = ds['t2m'].mean(dim=['valid_time']) - 273.15   # Convert to Celsius
    # Relative humidity calculation: Magnus-Tetens
    rh_mean = (100 * (np.exp(17.62 * ds['d2m'] / (243.12 + ds['d2m'])) / np.exp(17.62 * ds['t2m'] / (243.12 + ds['t2m'])))).mean(dim=['valid_time'])

    # Calculate sum of tp
    tp_sum = ds['tp'].sum(dim=['valid_time'])

    last_date = ds.valid_time.max().values.astype('datetime64[s]').astype(datetime).date()

    ds_agg = xr.Dataset(
        {
            'd2m_mean': d2m_mean,
            't2m_mean': t2m_mean,
            'rh_mean': rh_mean,
            'tp_sum': tp_sum
        },
        coords={
            'latitude': ds.latitude,
            'longitude': ds.longitude,
            'date': np.array([np.datetime64(last_date)])
        }
    )

    ds_agg['d2m_mean'].attrs['title'] = 'Avg Dew Point Temperature'
    ds_agg['d2m_mean'].attrs['units'] = 'degrees Celsius'
    ds_agg['t2m_mean'].attrs['title'] = 'Avg Temperature 2m'
    ds_agg['t2m_mean'].attrs['units'] = 'degrees Celsius'
    ds_agg['rh_mean'].attrs['title'] = 'Avg Relative Humidity'
    ds_agg['rh_mean'].attrs['units'] = '%'
    ds_agg['tp_sum'].attrs['title'] = 'Total precipitation'
    ds_agg['tp_sum'].attrs['units'] = 'm'  # for total precipitation

    ds_agg.to_netcdf(output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Aggregate ERA5 GRIB data daily.')
    parser.add_argument('--input_file', help='Input GRIB file path')
    parser.add_argument('--output_file', help='Output NetCDF file path')

    args = parser.parse_args()

    main(args.input_file, args.output_file)