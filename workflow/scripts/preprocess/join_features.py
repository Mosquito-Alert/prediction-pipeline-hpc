import argparse
import pandas as pd

try:
    sampling_effort_file_default = snakemake.input["sampling_effort"]
    bite_reports_file_default = snakemake.input["bite_reports"]
    era5_file_default = snakemake.input["era5"]
    landcover_file_default = snakemake.input["landcover"]
    output_file_default = snakemake.output[0]
except NameError:
    sampling_effort_file_default = None
    bite_reports_file_default = None
    era5_file_default = None
    landcover_file_default = None
    output_file_default = None


def main(sampling_effort_file: str, bite_reports_file: str, era5_file: str, landcover_file: str, output_file: str):
    df_se = pd.read_csv(sampling_effort_file, index_col=['h3_index', 'date'])
    df_bites = pd.read_csv(bite_reports_file, index_col=['h3_index', 'date'])
    df_era5 = pd.read_csv(era5_file, index_col=['h3_index', 'date'])
    df_landcover = pd.read_csv(landcover_file, index_col='h3_index')

    merged_df = df_se[df_se['SE'] > 0].merge(
        df_era5,
        left_index=True,
        right_index=True
    ).merge(
        df_bites,
        left_index=True,
        right_index=True,
        how='left'
    ).merge(
        df_landcover,
        left_index=True,
        right_index=True,
        how='left'
    ).fillna({
        'n_bite_reports': 0,
        'perc_discont_urban_fabric': 0,
        'perc_other_artificial': 0,
        'perc_agricultural': 0,
        'perc_other': 1
    }).reset_index()

    # Save the result to the output file
    merged_df.to_csv(output_file, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge CSV files on the 'h3_index' column.")
    parser.add_argument("--sampling_effort_file", type=str, required=(sampling_effort_file_default is None), default=sampling_effort_file_default)
    parser.add_argument("--bite_reports_file", type=str, required=(bite_reports_file_default is None), default=bite_reports_file_default)
    parser.add_argument("--era5_file", type=str, required=(era5_file_default is None), default=era5_file_default)
    parser.add_argument("--landcover_file", type=str, required=(landcover_file_default is None), default=landcover_file_default)
    parser.add_argument("--output_file", type=str, required=(output_file_default is None), default=output_file_default)
    args = parser.parse_args()

    main(
        sampling_effort_file=args.sampling_effort_file,
        bite_reports_file=args.bite_reports_file,
        era5_file=args.era5_file,
        landcover_file=args.landcover_file,
        output_file=args.output_file
    )