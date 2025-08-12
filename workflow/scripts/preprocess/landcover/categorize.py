import argparse
import pandas as pd

try:
    input_file_default = snakemake.input[0]
    output_file_default = snakemake.output[0]
except NameError:
    input_file_default = None
    output_file_default = None

def categorize(input_file: str, output_file: str):

    df = pd.read_csv(input_file)

    group_map = {
        "discont_urban_fabric": ['2'],
        "other_artificial": ['3', '5', '6', '7', '8', '9'],
        "agricultural": [str(i) for i in range(12, 23)],
    }
    all_cols = set(df.columns) - {'h3_index'}
    used_cols = set(sum(group_map.values(), []))
    group_map["other"] = list(all_cols - used_cols)

    # now create a new DataFrame with grouped sums
    grouped_df = pd.DataFrame()
    grouped_df["h3_index"] = df["h3_index"]

    for group, cols in group_map.items():
        grouped_df[group] = df[cols].sum(axis=1)

    grouped_df.to_csv(output_file, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Classify selected CLC categories from a GeoTIFF file.")
    parser.add_argument("-i", "--input_file", required=(input_file_default is None), default=input_file_default, help="Path to input .tif file")
    parser.add_argument("-o", "--output_file", required=(output_file_default is None), default=output_file_default, help="Path to output .tif file")
    args = parser.parse_args()

    categorize(input_file=args.input_file, output_file=args.output_file)