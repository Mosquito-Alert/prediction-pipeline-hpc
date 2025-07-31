include: "generate_bite_reports.smk"
include: "generate_era5_features.smk"
include: "generate_landcover.smk"
include: "generate_sampling_effort_features.smk"

rule generate_features:
    input:
        "outputs/features/bite_reports/{year}-{month}-{day}.csv",
        "outputs/features/era5/{year}-{month}-{day}.csv",
        "outputs/features/landcover.csv",
        "outputs/features/sampling_effort/{year}-{month}-{day}.csv"
    output:
        "outputs/features/all/{year}-{month}-{day}.csv"
    conda:
        "../envs/global.yaml"
    shell:
        "python3 workflow/scripts/preprocess/join_features.py --input_files {input} --output_file {output}"

from datetime import datetime, timedelta
def concat_all_features_inputs(wildcards):
    end_date = datetime(int(wildcards.year), int(wildcards.month), int(wildcards.day))
    start_date = datetime(2020, 1, 1)
    return [
        f"outputs/features/all/{(start_date + timedelta(days=i)).strftime('%Y-%m-%d')}.csv"
        for i in range((end_date - start_date).days + 1)
    ]

rule generate_train_dataset:
    input:
        input_files=concat_all_features_inputs
    output:
        "outputs/features/train_{year}-{month}-{day}.csv"
    conda:
        "../envs/global.yaml"
    shell:
        "python3 workflow/scripts/preprocess/generate_train_dataset.py --input_files {input.input_files} --output_file {output}"