rule download_era5:
    output:
        temp("data/era5/raw/{year}-{month}-{day}.zip", group_jobs=True)
    resources:
        runtime="3h"
    conda:
        "../envs/global.yaml"
    group:
        "era5_preprocess_{year}-{month}-{day}"
    params:
        date=lambda wildcards: f"{wildcards.year}-{wildcards.month}-{wildcards.day}",
        min_lon=config['era5']['min_lon'],
        min_lat=config['era5']['min_lat'],
        max_lon=config['era5']['max_lon'],
        max_lat=config['era5']['max_lat'],
    script:
        "../scripts/preprocess/era5/download.py"

rule convert_to_h3_era5:
    input:
        "data/era5/raw/{year}-{month}-{day}.zip"
    output:
        temp("data/era5/raw_h3/{year}-{month}-{day}.nc")
    conda:
        "../envs/global.yaml"
    group:
        "era5_preprocess_{year}-{month}-{day}"
    params:
        h3_res=config['h3_res'],
    script:
        "../scripts/preprocess/era5/convert_to_h3.py"

rule daily_stats_era5:
    input:
        "data/era5/raw_h3/{year}-{month}-{day}.nc"
    output:
        "data/era5/daily_stats/{year}-{month}-{day}.nc"
    conda:
        "../envs/global.yaml"
    group:
        "era5_preprocess_{year}-{month}-{day}"
    script:
        "../scripts/preprocess/era5/daily_stats.py"


from datetime import datetime, timedelta

def get_era5_data_preparation_inputs(wildcards):
    end_date = datetime(int(wildcards.year), int(wildcards.month), int(wildcards.day))
    start_date = end_date - timedelta(days=21)
    return [
        f"data/era5/daily_stats/{(start_date + timedelta(days=i)).strftime('%Y-%m-%d')}.nc"
        for i in range((end_date - start_date).days + 1)
    ]

rule generate_era5_features:
    input:
        get_era5_data_preparation_inputs
    output:
        "outputs/features/era5/{year}-{month}-{day}.csv"
    conda:
        "../envs/global.yaml"
    script:
        "../scripts/preprocess/era5/generate_features.py"