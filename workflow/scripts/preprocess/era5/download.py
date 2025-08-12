import argparse
import cdsapi
from datetime import datetime
from typing import Optional

CLIPPING_BUFFER = 0.5  # Buffer in degrees for clipping

# Try to import snakemake variables, if running within Snakemake
try:
    date_default = snakemake.params["date"]
    min_lon = snakemake.params['min_lon']
    min_lat = snakemake.params['min_lat']
    max_lon = snakemake.params['max_lon']
    max_lat = snakemake.params['max_lat']
    output_file_default = snakemake.output[0]
except NameError:
    date_default = None
    min_lon = None
    min_lat = None
    max_lon = None
    max_lat = None
    output_file_default = None

def main(date: datetime, output_file: str, min_lon: Optional[float] = None, min_lat: Optional[float] = None,
        max_lon: Optional[float] = None, max_lat: Optional[float] = None):
    c = cdsapi.Client()

    extra_requests = {}
    if None not in (min_lon, min_lat, max_lon, max_lat):
        extra_requests['area'] = [
            max_lat + CLIPPING_BUFFER, min_lon - CLIPPING_BUFFER,
            min_lat - CLIPPING_BUFFER, max_lon + CLIPPING_BUFFER
        ]
    elif any(v is not None for v in (min_lon, min_lat, max_lon, max_lat)):
        raise ValueError("If specifying bounds, all of min_lon, min_lat, max_lon, and max_lat must be provided.")

    c.retrieve(
        name='reanalysis-era5-single-levels',
        request={
            "product_type": ["reanalysis"],
            "variable": [
                "2m_dewpoint_temperature",
                "2m_temperature",
                "total_precipitation"
            ],
            "year": str(date.year),
            "month": str(date.month).zfill(2),
            "day": [str(date.day).zfill(2)],
            "time": [
                "00:00", "01:00", "02:00",
                "03:00", "04:00", "05:00",
                "06:00", "07:00", "08:00",
                "09:00", "10:00", "11:00",
                "12:00", "13:00", "14:00",
                "15:00", "16:00", "17:00",
                "18:00", "19:00", "20:00",
                "21:00", "22:00", "23:00"
            ],
            "data_format": "netcdf",
            "download_format": "unarchived",
            **extra_requests
        },
        target=output_file,
    )

if __name__ == "__main__":
    def valid_date(s):
        try:
            return datetime.strptime(s, "%Y-%m-%d").date()
        except ValueError:
            raise argparse.ArgumentTypeError(f"Not a valid date: '{s}'. Expected format: YYYY-MM-DD.")

    parser = argparse.ArgumentParser(description="Download ERA5 Land data via CDS API.")
    parser.add_argument('--date', type=valid_date, required=(date_default is None), default=date_default, help="Date in YYYY-MM-DD format")
    parser.add_argument('--output_file', type=str, required=(output_file_default is None), default=output_file_default, help='Output filename')
    parser.add_argument('--min_lon', type=float, default=min_lon, help='Minimum longitude')
    parser.add_argument('--min_lat', type=float, default=min_lat, help='Minimum latitude')
    parser.add_argument('--max_lon', type=float, default=max_lon, help='Maximum longitude')
    parser.add_argument('--max_lat', type=float, default=max_lat, help='Maximum latitude')

    args = parser.parse_args()

    main(
        date=args.date,
        output_file=args.output_file,
        min_lon=args.min_lon,
        min_lat=args.min_lat,
        max_lon=args.max_lon,
        max_lat=args.max_lat
    )