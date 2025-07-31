import argparse
from datetime import datetime, timezone
from flatten_json import flatten
import pandas as pd

import mosquito_alert

def main(date: str, output_file: str):
    bites = []

    date = datetime.strptime(date, "%Y-%m-%d")

    with mosquito_alert.ApiClient() as api_client:
        bites_api = mosquito_alert.BitesApi(api_client)
        page = 1
        created_at_after = datetime.combine(date, datetime.min.time(), tzinfo=timezone.utc)
        created_at_before = datetime.combine(date, datetime.max.time(), tzinfo=timezone.utc)
        while True:
            api_response = bites_api.list(
                page=page,
                page_size=100,
                created_at_after=created_at_after,
                created_at_before=created_at_before,
                order_by=['created_at'],
            )
            bites.extend(api_response.results)
            if not api_response.next:
                break
            page += 1
    df = pd.DataFrame([flatten(b.model_dump(by_alias=True)) for b in bites])
    df.to_csv(output_file, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download mosquito bites for a specific date.")
    parser.add_argument("--date", type=str, help="Date in YYYY-MM-DD format")
    parser.add_argument("--output_file", type=str, help="Output CSV file path")
    args = parser.parse_args()

    main(args.date, args.output_file)