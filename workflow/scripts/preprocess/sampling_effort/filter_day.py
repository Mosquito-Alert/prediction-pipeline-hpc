import argparse
import pandas as pd

# SAMPLING_EFFORT_URL = 'https://github.com/Mosquito-Alert/sampling_effort_data/raw/main/sampling_effort_daily_cellres_025.csv.gz'

def main(input_file: str, output_file: str, date: str):

    df = pd.read_csv(
        input_file,
        date_format='%Y-%m-%d',
        parse_dates=['date']
    )

    # Filter by date
    df = df[df['date'] == pd.to_datetime(date)]

    df.drop(columns=['TigacellID'], inplace=True)
    df.rename(columns={'masked_lon': 'longitude', 'masked_lat': 'latitude'}, inplace=True)
    
    df['n_reports_total'] = df.filter(like='n_reports_').sum(axis=1)
    df['n_reporters_total'] = df.filter(like='n_reporters_').sum(axis=1)
    
    df.to_csv(output_file, index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process sampling effort data.')
    parser.add_argument('--input_file', type=str, required=True, help='Input CSV file path')
    parser.add_argument('--output_file', type=str, required=True, help='Output CSV file path')
    parser.add_argument('--date', type=str, required=False, help='Filter data by date (YYYY-MM-DD)')

    args = parser.parse_args()
    main(args.input_file, args.output_file, args.date)