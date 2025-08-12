import argparse
import pandas as pd

try:
    input_file_default = snakemake.input[0]
    output_file_default = snakemake.output[0]
except NameError:
    input_file_default = None
    output_file_default = None

def generate_features(input_file: str, output_file: str):

    df = pd.read_csv(input_file, index_col='h3_index')

    columns = df.columns

    df["total_count"] = df.sum(axis=1)

    for col in columns:
        df[f"perc_{col}"] = df[col] / df['total_count']

    df.to_csv(output_file, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Classify selected CLC categories from a GeoTIFF file.")
    parser.add_argument("-i", "--input_file", required=(input_file_default is None), default=input_file_default, help="Path to input .tif file")
    parser.add_argument("-o", "--output_file", required=(output_file_default is None), default=output_file_default, help="Path to output .tif file")
    args = parser.parse_args()

    generate_features(input_file=args.input_file, output_file=args.output_file)