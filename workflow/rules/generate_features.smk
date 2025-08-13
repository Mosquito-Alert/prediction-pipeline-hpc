include: "generate_bite_reports.smk"
include: "generate_era5_features.smk"
include: "generate_landcover.smk"
include: "generate_sampling_effort_features.smk"

rule generate_features:
    input:
        sampling_effort="outputs/features/sampling_effort/{year}-{month}-{day}.csv",
        bite_reports="outputs/features/bite_reports/{year}-{month}-{day}.csv",
        era5="outputs/features/era5/{year}-{month}-{day}.csv",
        landcover="outputs/features/landcover/h3_landcover_stats.csv",
    output:
        "outputs/features/merged/{year}-{month}-{day}.csv"
    conda:
        "../envs/global.yaml"
    script:
        "../scripts/preprocess/join_features.py"

from datetime import datetime, timedelta
def concat_all_features_inputs(wildcards):
    end_date = datetime(int(wildcards.year), int(wildcards.month), int(wildcards.day))
    start_date = datetime(2020, 1, 1)
    return [
        f"outputs/features/merged/{(start_date + timedelta(days=i)).strftime('%Y-%m-%d')}.csv"
        for i in range((end_date - start_date).days + 1)
    ]

rule generate_train_dataset:
    input:
        input_files=concat_all_features_inputs
    output:
        "outputs/features/train/{year}-{month}-{day}.csv"
    conda:
        "../envs/global.yaml"
    script:
        "../scripts/preprocess/generate_train_dataset.py"