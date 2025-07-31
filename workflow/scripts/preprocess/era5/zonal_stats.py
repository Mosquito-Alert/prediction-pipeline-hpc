import argparse
import geopandas as gpd
import numpy as np
import regionmask
import xarray as xr

def main(input_file, vector_file, output_file):
    ds = xr.open_dataset(input_file)
    ds = ds.interp(
        latitude=np.linspace(ds.latitude.min().item(), ds.latitude.max().item(), len(ds.latitude)),
        longitude=np.linspace(ds.longitude.min().item(), ds.longitude.max().item(), len(ds.longitude))
    )

    gdf = gpd.read_file(vector_file)

    # Alternatives:
    # See: https://github.com/regionmask/regionmask/issues/225#issuecomment-915033614
    # See: https://github.com/regionmask/regionmask/issues/304#issuecomment-1007938765
    # See: https://xarray-spatial.readthedocs.io/en/latest/user_guide/zonal.html
    # See: https://medium.com/@lubomirfranko/perform-zonal-statistics-on-climate-data-with-python-4df6b7e5a171
    # See: https://pangeo-xesmf.readthedocs.io/en/latest/notebooks/Spatial_Averaging.html

    regions = regionmask.from_geopandas(gdf, names='code')

    mask = regions.mask_3D_frac_approx(ds.longitude, ds.latitude, drop=False)
    region_means_ds = ds.where(mask > 0).weighted(mask).mean(dim=("latitude", "longitude"))
    # Add date dimension in the result
    region_means_ds = region_means_ds.expand_dims(date=ds.date)
    region_means_ds.to_netcdf(output_file)

    # region_means_ds = ds.where(mask > 0).groupby('region').weighted(mask).mean(["latitude", "longitude"])

    # # mask = regionmask.mask_3D_geopandas(gdf, ds.longitude, ds.latitude, drop=False)

    
    # mask_all = regionmask.core.mask._mask_rasterize(ds.longitude, ds.latitude, regions.polygons, regions.numbers, all_touched=True, as_3D=True)
    # regionmask.core.mask._mask_to_dataarray(mask_all, ds.longitude, ds.latitude)

    # region_means_ds = ds.where(mask).groupby('region').weighted(mask).mean(["latitude","longitude"])

    # df = region_means_ds.to_dataframe()

    # gdf.join(df).to_csv('region_means.csv')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute region means from ERA5 dataset using vector mask.")
    parser.add_argument("--input_file", type=str, required=True, help="Input NetCDF file path")
    parser.add_argument("--vector_file", type=str, required=True, help="Vector file path (e.g. GeoPackage)")
    parser.add_argument("--output_file", type=str, required=True, help="Output NetCDF file path")

    args = parser.parse_args()
    main(args.input_file, args.vector_file, args.output_file)