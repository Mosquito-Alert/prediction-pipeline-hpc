import argparse
import pandas as pd
from typing import List

try:
    input_files_default = snakemake.input
    output_file_default = snakemake.output[0]
except NameError:
    input_files_default = None
    output_file_default = None

def main(input_files: List[str], output_file: str):
    # Read the input files
    dfs = [pd.read_csv(file) for file in input_files]

    # Concatenate all dataframes vertically
    concat_df = pd.concat(dfs, ignore_index=True)

    # Filter out rows where SE is 0 or NA
    filtered_df = concat_df[(concat_df['SE'].notna()) & (concat_df['SE'] != 0)]

    counts = filtered_df['h3_index'].value_counts()
    valid_indices = counts[counts >= 3].index
    df_filtered = filtered_df[filtered_df['h3_index'].isin(valid_indices)]

    # Save the result to the output file
    df_filtered.to_csv(output_file, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Concat CSV feature files.")
    parser.add_argument("--input_files", nargs='+', required=(input_files_default is None), default=input_files_default, help="List of input CSV files (space-separated).")
    parser.add_argument("--output_file", type=str, required=(output_file_default is None), default=output_file_default, help="Output CSV file")

    args = parser.parse_args()

    main(args.input_files, args.output_file)