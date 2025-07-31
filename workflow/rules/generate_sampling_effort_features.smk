rule download_sampling_effort:
    output:
        temp("data/sampling_effort/sampling_effort_daily_cellres_025.csv.gz")
    resources:
        api_calls=5
    params:
        url=config['get_sampling_effort']['url']
    shell:
        "curl -L -o {output} {params.url}"

rule daily_sampling_effort:
    input:
        "data/sampling_effort/sampling_effort_daily_cellres_025.csv.gz"
    output:
        "data/sampling_effort/{year}-{month}-{day}.csv"
    conda:
        "../envs/global.yaml"
    shell:
        "python3 workflow/scripts/preprocess/sampling_effort/filter_day.py --input_file {input} --date {wildcards.year}-{wildcards.month}-{wildcards.day} --output_file {output}"

rule zonal_stats_sampling_effort:
    input:
        sampling_effort="data/sampling_effort/{year}-{month}-{day}.csv",
        vector="data/gadm_410_esp_simplified.gpkg"
    output:
        "outputs/features/sampling_effort/{year}-{month}-{day}.csv"
    conda:
        "../envs/global.yaml"
    shell:
        "python3 workflow/scripts/preprocess/sampling_effort/zonal_stats.py --input_file {input.sampling_effort} --vector_file {input.vector} --output_file {output}"