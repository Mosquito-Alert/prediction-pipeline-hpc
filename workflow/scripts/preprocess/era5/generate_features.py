import argparse
import h3
import xarray as xr

try:
    input_files_default = snakemake.input
    output_file_default = snakemake.output[0]
except NameError:
    input_files_default = None
    output_file_default = None

def main(input_files, output_file):
    # Open multiple NetCDF files and concatenate them along a new dimension
    ds = xr.open_mfdataset(input_files, combine='nested', concat_dim="date", parallel=True)

    min_t2m_21d = ds['t2m_mean'].min(dim='date')
    last_date_ds = ds.sel(date=ds.date.max())
    last_date_ds = last_date_ds.assign(min_t2m_21d=min_t2m_21d)

    df = last_date_ds.to_dataframe().reset_index()
    df['h3_index'] = df['cell_ids'].apply(h3.int_to_str)
    df = df.drop(['number', 'cell_ids'], axis=1)

    # Save the combined dataset to a new NetCDF file
    df.to_csv(output_file, index=False)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute region means from ERA5 dataset using vector mask.")
    parser.add_argument("--input_files", nargs='+', required=(input_files_default is None), default=input_files_default, help="Input NetCDF file to combine stats from")
    parser.add_argument("--output_file", type=str, required=(output_file_default is None), default=output_file_default, help="Output CSV file")

    args = parser.parse_args()
    main(args.input_files, args.output_file)