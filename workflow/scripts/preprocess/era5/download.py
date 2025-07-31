import argparse
import fiona
import cdsapi

CLIPPING_BUFFER = 0.1  # Buffer in degrees for clipping

def main():
    parser = argparse.ArgumentParser(description="Download ERA5 Land data via CDS API.")
    parser.add_argument('--year', type=int, required=True, help='Year')
    parser.add_argument('--month', type=int, required=True, help='Month')
    parser.add_argument('--day', type=int, required=True, help='Day')
    parser.add_argument('--output_file', type=str, required=True, help='Output filename')
    parser.add_argument('--clip_vector', type=str, required=False, help='Clip vector')

    args = parser.parse_args()

    c = cdsapi.Client()

    extra_requests = {}
    if args.clip_vector:
        with fiona.open(args.clip_vector) as vector_file:
            minx, miny, maxx, maxy = vector_file.bounds
            extra_requests['area'] = [
                maxy + CLIPPING_BUFFER, minx - CLIPPING_BUFFER,
                miny - CLIPPING_BUFFER, maxx + CLIPPING_BUFFER
            ]

    c.retrieve(
        name='reanalysis-era5-land',
        request={
            "variable": [
                "2m_dewpoint_temperature",
                "2m_temperature",
                "total_precipitation"
            ],
            "year": str(args.year),
            "month": str(args.month).zfill(2),
            "day": [str(args.day).zfill(2)],
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
        target=args.output_file,
    )

if __name__ == "__main__":
    main()