import argparse
from datetime import datetime
import numpy as np
import os
import tempfile
import xarray as xr
import zipfile

try:
    input_file_default = snakemake.input[0]
    output_file_default = snakemake.output[0]
except NameError:
    input_file_default = None
    output_file_default = None

# TODO: check why there's more than 1 day.

def main(input_file: str, output_file: str):

    with tempfile.TemporaryDirectory() as tmpdir:
        # Unzip the input file
        with zipfile.ZipFile(input_file, 'r') as zip_ref:
            zip_ref.extractall(tmpdir)

        # Find the unzipped NetCDF file
        extracted_files = [os.path.join(tmpdir, f) for f in os.listdir(tmpdir) if f.endswith('.nc')]
        if not extracted_files:
            raise FileNotFoundError("No NetCDF files found in the unzipped content.")

        ds = xr.open_mfdataset(extracted_files)

    # Calculate mean of d2m and t2m
    d2m_celsius = ds['d2m'] - 273.15  # Convert from Kelvin to Celsius
    t2m_celsius = ds['t2m'] - 273.15  # Convert from Kelvin to Celsius
    # Relative humidity calculation: Magnus-Tetens
    # rh_mean = (100 * (np.exp(17.62 * d2m_celsius / (243.12 + d2m_celsius)) / np.exp(17.62 * t2m_celsius / (243.12 + t2m_celsius)))).mean(dim=['valid_time'])
    # See: https://earthscience.stackexchange.com/questions/16570/how-to-calculate-relative-humidity-from-temperature-dew-point-and-pressure
    m = 7.591386
    Tn = 240.7263
    rh_mean = 100 * 10 ** (m * ((d2m_celsius/(d2m_celsius+Tn))-(t2m_celsius/(t2m_celsius+Tn)))).mean(dim=['valid_time'])

    d2m_mean = d2m_celsius.mean(dim=['valid_time'])
    t2m_mean = t2m_celsius.mean(dim=['valid_time'])
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
    parser.add_argument('--input_file', type=str, required=(input_file_default is None), default=input_file_default, help='Input GRIB file path')
    parser.add_argument('--output_file', type=str, required=(output_file_default is None), default=output_file_default, help='Output NetCDF file path')

    args = parser.parse_args()

    main(args.input_file, args.output_file)