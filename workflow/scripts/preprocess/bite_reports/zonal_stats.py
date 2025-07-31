import argparse
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point

def main(vector_file: str, input_file: str, output_file: str):
    gdf = gpd.read_file(vector_file)
    gdf.drop(columns=['name'], inplace=True)
    gdf.rename(columns={'code': 'id'}, inplace=True)

    df = pd.read_csv(input_file)
    df = df[df['counts_total'] > 0]
    geometry = [Point(xy) for xy in zip(df['location_point_longitude'], df['location_point_latitude'])]
    gdf_bites = gpd.GeoDataFrame(df, geometry=geometry, crs='EPSG:4326')
    gdf_bites.drop(columns=['location_point_longitude', 'location_point_latitude'], inplace=True)
    del df

    # Spatial join and aggregation
    df_count = gpd.sjoin(gdf_bites, gdf, how="inner", predicate="within").groupby("id", as_index=False).size()
    df_count.rename(columns={'size': 'n_bite_reports'}, inplace=True)

    # Merge and export
    gdf_result = gdf.merge(df_count, left_on='id', right_on='id', how='left')
    gdf_result['n_bite_reports'] = gdf_result['n_bite_reports'].fillna(0)
    gdf_result['has_bite_reports'] = gdf_result['n_bite_reports'] > 0

    gdf_result.drop(columns=['geometry'], inplace=True)
    gdf_result.to_csv(output_file, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Aggregate bites reports spatially.")
    parser.add_argument("--input_file", help="Path to input CSV file with bites reports")
    parser.add_argument("--vector_file", help="Path to vector file (e.g., GPKG or shapefile)")
    parser.add_argument("--output_file", help="Path to output CSV file")

    args = parser.parse_args()
    main(args.vector_file, args.input_file, args.output_file)