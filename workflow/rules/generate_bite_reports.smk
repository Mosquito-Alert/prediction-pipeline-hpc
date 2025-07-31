rule download_bite_reports:
    output:
        "data/bite_reports/{year}-{month}-{day}.csv"
    retries: 3
    resources:
        api_calls=5
    conda:
        "../envs/global.yaml"
    shell:
        "python3 workflow/scripts/preprocess/bite_reports/download.py --date {wildcards.year}-{wildcards.month}-{wildcards.day} --output_file {output}"

rule zonal_stats_bite_reports:
    input:
        bite_reports="data/bite_reports/{year}-{month}-{day}.csv",
        vector="data/gadm_410_esp_simplified.gpkg"
    output:
        "outputs/features/bite_reports/{year}-{month}-{day}.csv"
    conda:
        "../envs/global.yaml"
    shell:
        "python3 workflow/scripts/preprocess/bite_reports/zonal_stats.py --input_file {input.bite_reports} --vector_file {input.vector} --output_file {output}"