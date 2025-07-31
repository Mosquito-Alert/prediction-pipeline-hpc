rule categorize_landcover:
    input:
        "data/corine_landcover/U2018_CLC2018_V2020_20u1.tif"
    output:
        "data/corine_landcover/landcover_categorized.tif"
    shell:
        "python3 workflow/scripts/preprocess/landcover/categorize.py --input_file {input} --output_file {output}"

rule landcover_zonal_stats:
    input:
        landcover="data/corine_landcover/landcover_categorized.tif",
        vector="data/gadm_410_esp_simplified.gpkg"
    output:
        "outputs/features/landcover.csv"
    shell:
        "python3 workflow/scripts/preprocess/landcover/zonal_stats.py --vector_file {input.vector} --raster_file {input.landcover} --output_csv {output}"