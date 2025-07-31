import argparse
import pandas as pd
from typing import List

def main(input_files: List[str], output_file: str):
    # Read the input files
    dfs = [pd.read_csv(file) for file in input_files]

    # Merge all dataframes on the 'id' column
    merged_df = dfs[0]
    for df in dfs[1:]:
        merged_df = pd.merge(merged_df, df, on='id', how='inner')

    # Save the result to the output file
    merged_df.to_csv(output_file, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge CSV files on the 'id' column.")
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