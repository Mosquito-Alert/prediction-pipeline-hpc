import argparse
import h3
import pandas as pd

try:
    input_file_default = snakemake.input[0]
    h3_resolution = snakemake.params["h3_res"]
    output_file_default = snakemake.output[0]
except NameError:
    input_file_default = None
    h3_resolution = None
    output_file_default = None

def main(input_file: str, output_file: str, h3_res: int):
    df = pd.read_csv(input_file)
    df = df[df['counts_total'] > 0]

    if not df.empty:
        df['h3_index'] = df.apply(
            lambda row: h3.latlng_to_cell(
                row['location_point_latitude'],
                row['location_point_longitude'],
                h3_res
            ),
            axis=1
        )
        df_agg = df.groupby('h3_index').size().reset_index(name='n_bite_reports')
    else:
        df_agg = pd.DataFrame(columns=['h3_index', 'n_bite_reports'])

    df_agg.to_csv(output_file, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Aggregate bites reports spatially.")
    parser.add_argument("--input_file", type=str, required=(input_file_default is None), default=input_file_default, help="Path to input CSV file with bites reports")
    parser.add_argument('--h3_resolution', type=int, required=(h3_resolution is None), default=h3_resolution, help='The H3 resolution to use')
    parser.add_argument("--output_file", type=str, required=(output_file_default is None), default=output_file_default, help="Path to output CSV file")

    args = parser.parse_args()
    main(args.input_file, args.output_file, args.h3_resolution)