import argparse
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point

def main(vector_file: str, input_file: str, output_file: str):
    gdf = gpd.read_file(vector_file)
    gdf.drop(columns=['name'], inplace=True)
    gdf.rename(columns={'code': 'id'}, inplace=True)

    df = pd.read_csv(input_file)
    geometry = [Point(xy) for xy in zip(df['longitude'], df['latitude'])]
    gdf_sampling_effort = gpd.GeoDataFrame(df, geometry=geometry, crs='EPSG:4326')
    gdf_sampling_effort.drop(columns=['longitude', 'latitude'], inplace=True)
    del df

    # Spatial join and aggregation
    df_grouped = gpd.sjoin(gdf_sampling_effort, gdf, how="inner", predicate="within").groupby("id", as_index=False)
    df_agg = df_grouped.agg(
        SE=('SE', lambda x: 1 - np.prod(1 - x)),
        n_reports_total=('n_reports_total', 'sum'),
        n_reporters_total=('n_reporters_total', 'sum')
    )

    # Merge and export
    gdf_result = gdf.merge(df_agg, left_on='id', right_on='id', how='left')
    gdf_result.drop(columns=['geometry'], inplace=True)
    gdf_result.to_csv(output_file, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Aggregate sampling effort spatially.")
    parser.add_argument("--input_file", help="Path to input CSV file with sampling effort")
    parser.add_argument("--vector_file", help="Path to vector file (e.g., GPKG or shapefile)")
    parser.add_argument("--output_file", help="Path to output CSV file")

    args = parser.parse_args()
    main(args.vector_file, args.input_file, args.output_file)