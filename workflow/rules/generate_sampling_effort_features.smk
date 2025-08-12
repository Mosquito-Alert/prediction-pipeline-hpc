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
    params:
        date=lambda wildcards: f"{wildcards.year}-{wildcards.month}-{wildcards.day}",
    script:
        "workflow/scripts/preprocess/sampling_effort/filter_date.py"

rule convert_to_h3_sampling_effort:
    input:
        "data/sampling_effort/{year}-{month}-{day}.csv"
    output:
        "outputs/features/sampling_effort/{year}-{month}-{day}.csv"
    conda:
        "../envs/global.yaml"
    params:
        h3_res=config['h3_res'],
    script:
        "workflow/scripts/preprocess/sampling_effort/convert_to_h3.py"
