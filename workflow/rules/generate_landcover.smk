rule convert_to_h3_landcover:
    input:
        "data/corine_landcover/U2018_CLC2018_V2020_20u1.tif"
    output:
        "data/corine_landcover/h3_landcover.csv"
    conda:
        "../envs/global.yaml"
    params:
        h3_res=config['h3_res'],
    script:
        "../scripts/preprocess/landcover/convert_to_h3.py"

rule categorize_landcover:
    input:
        "data/corine_landcover/h3_landcover.csv"
    output:
        "data/corine_landcover/h3_landcover_categorized.csv"
    conda:
        "../envs/global.yaml"
    script:
        "../scripts/preprocess/landcover/categorize.py"

rule generate_features_landcover:
    input:
        "data/corine_landcover/h3_landcover_categorized.csv"
    output:
        "outputs/features/landcover/h3_landcover_stats.csv"
    conda:
        "../envs/global.yaml"
    script:
        "../scripts/preprocess/landcover/generate_features.py"