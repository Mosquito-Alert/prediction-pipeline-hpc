import argparse
from collections import defaultdict
import h3
import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import transform
from rasterio.windows import Window

try:
    input_file_default = snakemake.input[0]
    h3_resolution = snakemake.params["h3_res"]
    output_file_default = snakemake.output[0]
except NameError:
    input_file_default = None
    h3_resolution = None
    output_file_default = None

BLOCK_HEIGHT = 1024
BLOCK_WIDTH = 1024

def main(input_file: str, output_file: str, h3_res: int):
    h3_counts = defaultdict(lambda: defaultdict(int))

    with rasterio.open(input_file) as src:
        for row_off in range(0, src.height, BLOCK_HEIGHT):
            for col_off in range(0, src.width, BLOCK_WIDTH):
                win = Window(
                    col_off,
                    row_off,
                    min(BLOCK_WIDTH, src.width - col_off),
                    min(BLOCK_HEIGHT, src.height - row_off)
                )
                data = src.read(1, window=win, masked=True)
                if data.mask.all():
                    continue
                # Flatten indices of valid pixels
                valid_pixels = ~data.mask
                rows, cols = np.where(valid_pixels)
                values = data[rows, cols]
                # Map pixel indices to global raster coords
                global_rows = rows + row_off
                global_cols = cols + col_off
                # Get lat/lon of pixels
                xs, ys = src.xy(global_rows, global_cols)  # xs=lon, ys=lat
                lons, lats = transform(src.crs, 'EPSG:4326', xs, ys)
                h3_indexes = np.vectorize(lambda lat, lon: h3.latlng_to_cell(lat, lon, h3_res))(lats, lons)
                for h, v in zip(h3_indexes, values):
                    h3_counts[h][v] += 1

    df = pd.DataFrame.from_dict(h3_counts, orient='index')
    df = df.fillna(0)
    df = df.reset_index()
    df = df.rename(columns={'index':'h3_index'})
    df.to_csv(output_file, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Aggregate ERA5 GRIB data daily.')
    parser.add_argument('--input_file', type=str, required=(input_file_default is None), default=input_file_default, help='Input GRIB file path')
    parser.add_argument('--output_file', type=str, required=(output_file_default is None), default=output_file_default, help='Output NetCDF file path')
    parser.add_argument('--h3_resolution', type=int, required=(h3_resolution is None), default=h3_resolution, help='The H3 resolution to use')

    args = parser.parse_args()

    main(args.input_file, args.output_file, args.h3_resolution)