import argparse
import pandas as pd
from typing import List

def main(input_files: List[str], output_file: str):
    # Read the input files
    dfs = [pd.read_csv(file) for file in input_files]

    # Concatenate all dataframes vertically
    concat_df = pd.concat(dfs, ignore_index=True)

    # Filter out rows where SE is 0 or NA
    filtered_df = concat_df[(concat_df['SE'].notna()) & (concat_df['SE'] != 0)]

    # Save the result to the output file
    filtered_df.to_csv(output_file, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Concat CSV feature files.")
    parser.add_argument(
        "--input_files",
        nargs='+',
        type=str,
        required=True,
        help="List of input CSV files (space-separated)."
    )
    parser.add_argument(
        "--output_file",
        type=str,
        required=True,
        help="Path to output CSV file."
    )
    args = parser.parse_args()

    main(args.input_files, args.output_file)