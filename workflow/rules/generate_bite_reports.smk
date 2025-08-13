rule download_bite_reports:
    output:
        "data/bite_reports/{year}-{month}-{day}.csv"
    retries: 3
    conda:
        "../envs/global.yaml"
    group:
        "bite_preprocess_{year}-{month}-{day}"
    params:
        date=lambda wildcards: f"{wildcards.year}-{wildcards.month}-{wildcards.day}",
    script:
        "../scripts/preprocess/bite_reports/download.py"

rule convert_to_h3_bite_reports:
    input:
        "data/bite_reports/{year}-{month}-{day}.csv",
    output:
        "outputs/features/bite_reports/{year}-{month}-{day}.csv"
    conda:
        "../envs/global.yaml"
    group:
        "bite_preprocess_{year}-{month}-{day}"
    params:
        h3_res=config['h3_res'],
    script:
        "../scripts/preprocess/bite_reports/convert_to_h3.py"