rule download_era5:
    input:
        vector="data/gadm_410_esp_simplified.gpkg"
    output:
        temp("data/era5/{year}-{month}-{day}.grib", group_jobs=True)
    resources:
        api_calls=5,
        cpus_per_task=1
    conda:
        "../envs/global.yaml"
    params:
        date=lambda wildcards: f"{wildcards.year}-{wildcards.month}-{wildcards.day}",
    script:
        "../scripts/preprocess/era5/download.py"

rule daily_stats_era5:
    input:
        "data/era5/{year}-{month}-{day}.nc"
    output:
        "data/era5/daily_stats/{year}-{month}-{day}.nc"
    conda:
        "../envs/global.yaml"
    script:
        "../scripts/preprocess/era5/daily_stats.py"

rule zonal_stats_era5:
    input:
        era5="data/era5/daily_stats/{year}-{month}-{day}.nc",
        vector="data/gadm_410_esp_simplified.gpkg"
    output:
        "data/era5/zonal_stats/{year}-{month}-{day}.nc"
    conda:
        "../envs/global.yaml"
    resources:
        mem_mb=30000,
    script:
        "../scripts/preprocess/era5/zonal_stats.py"

from datetime import datetime, timedelta

def get_era5_data_preparation_inputs(wildcards):
    end_date = datetime(int(wildcards.year), int(wildcards.month), int(wildcards.day))
    start_date = end_date - timedelta(days=21)
    return [
        f"data/era5/zonal_stats/{(start_date + timedelta(days=i)).strftime('%Y-%m-%d')}.nc"
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