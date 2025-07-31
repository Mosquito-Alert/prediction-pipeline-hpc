import argparse
import geopandas as gpd
import pandas as pd
import rasterio
from rasterstats import gen_zonal_stats

CRS = 'EPSG:4326'

def compute_zonal_stats(vector_file: str, raster_file: str, output_csv: str, vector_layer: int = 0, raster_band: int = 1):
    raster = rasterio.open(raster_file, "r")

    gdf = gpd.read_file(vector_file, layer=vector_layer)
    gdf.drop(columns=['name'], inplace=True)

    # Generate zonal statistics
    zs = gen_zonal_stats(
        vectors=gdf.to_crs(raster.crs),
        raster=raster.read(raster_band),
        affine=raster.transform,
        nodata=raster.nodata,
        all_touched=True,
        categorical=True,
        category_map={
            0: 'other',
            1: 'discont_urban_fabric',
            2: 'other_artificial',
            3: 'agricultural'
        },
        prefix='count_'
    )
    zs_df = pd.DataFrame(zs)

    # Compute the percentage of each category
    count_cols = [col for col in zs_df.columns if col.startswith('count_')]
    zs_df[count_cols] = zs_df[count_cols].fillna(0)
    zs_df['total_count'] = zs_df[count_cols].sum(axis=1)

    for col in count_cols:
        perc_col = col.replace('count_', 'perc_')
        zs_df[perc_col] = (zs_df[col] / zs_df['total_count']) * 100 if zs_df['total_count'].any() else 0

    # Join with vector data
    gdf_with_stats = pd.concat([gdf.reset_index(drop=True), zs_df], axis=1)

    gdf_with_stats = gdf_with_stats.rename(columns={'code': 'id'})

    # Save to CSV
    gdf_with_stats.drop(columns="geometry").to_csv(output_csv, index=False)


def main():
    parser = argparse.ArgumentParser(description="Compute zonal statistics for vector polygons over a raster.")
    parser.add_argument('--vector_file', type=str, required=True, help="Path to the vector file (e.g., .gpkg)")
    parser.add_argument('--raster_file', type=str, required=True, help="Path to the raster file (e.g., .tif)")
    parser.add_argument('--output_csv', type=str, required=True, help="Output CSV file to save the results.")
    parser.add_argument('--vector_layer', type=int, default=0, help="Vector layer index (default: 0).")
    parser.add_argument('--raster_band', type=int, default=1, help="Raster band index (default: 1).")

    args = parser.parse_args()
    compute_zonal_stats(
        vector_file=args.vector_file,
        raster_file=args.raster_file,
        output_csv=args.output_csv,
        vector_layer=args.vector_layer,
        raster_band=args.raster_band
    )

if __name__ == "__main__":
    main()