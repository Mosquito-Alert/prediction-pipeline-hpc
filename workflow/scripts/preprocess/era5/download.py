import argparse
import cdsapi
from datetime import datetime
import fiona
from typing import Optional

CLIPPING_BUFFER = 0.5  # Buffer in degrees for clipping

# Try to import snakemake variables, if running within Snakemake
try:
    date_default = snakemake.params["date"]
    vector_file_deafult = snakemake.input["vector"]
    output_file_default = snakemake.output[0]
except NameError:
    date_default = None
    vector_file_deafult = None
    output_file_default = None

def main(date: datetime, output_file: str, clip_vector: Optional[str] = None):
    c = cdsapi.Client()

    extra_requests = {}
    if clip_vector:
        with fiona.open(clip_vector) as vector_file:
            minx, miny, maxx, maxy = vector_file.bounds
            extra_requests['area'] = [
                maxy + CLIPPING_BUFFER, minx - CLIPPING_BUFFER,
                miny - CLIPPING_BUFFER, maxx + CLIPPING_BUFFER
            ]

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
            "data_format": "grib",
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
    parser.add_argument('--clip_vector', type=str, required=(vector_file_deafult is None), default=vector_file_deafult, help='Clip vector')

    args = parser.parse_args()

    main(date=args.date, output_file=args.output_file, clip_vector=args.clip_vector)