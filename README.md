/home/ccerecedo/Bites_MA_Spain/

Nightly:
-> /leov1/bites_ma_spain/data/era5/YYYY-MM-DD.grib		(@epou)


R grib2csv.R
	--input_dir /leov1/bites_ma_spain/data/era5/
	--format YYYY-MM-DD.grid
	--boundary_mask /leov1/botes_ma_spain/data/boundaries/gadm41_ESP_simplified.gpkg
	--boundary_column_id code
	--date YYYY-MM-DD
	--output /leov1/bites_ma_spain/data/era5/YYYY-MM-DD.csv
	
	
	https://cran.r-project.org/web/packages/argparse/refman/argparse.html

		library(argparse)

		parser <- ArgumentParser()
		parser$add_argument("--input", help = "Path to data file", required = TRUE)
		parser$add_argument("--output", help = "Path to store resulting data file", required = TRUE)
		parser$add_argument("--boundary_mask", help = "Path to get vector boundaries to apply mask", required = TRUE)

		args <- parser$parse_args()

		print(args$input)
		print(args$output)
		print(args$boundary_mask)

	-> /leov1/bites_ma_spain/data/era5/YYYY-MM-DD.csv     (@catu, poner tambien los lags aqui)

	
Download Bites from API -> (@epou, mirar CSV que vaig enviar a la maria)
	
-> /leov1/bites_ma_spain/data/bites_ma/YYYY-MM-DD.csv

-> /leov1/bites_ma_spain/data/boundaries/gadm41_ESP_simplified.gpkg


/leov1/bites_ma_spain/data/corine_landcover/u2018.....tif
/leov1/bites_ma_spain/data/corine_landcover/u2018.....csv   (@catu, optional script que acepta argumentos de entrada)

(preparing_data_mun_template_bites)
R merge.R 
	--bites_csv /leov1/bites_ma_spain/data/bites_ma/YYYY-MM-DD.csv
	--boundary_mask /leov1/botes_ma_spain/data/boundaries/gadm41_ESP_simplified.gpkg 
	--boundary_column_id code 
	--landcover_csv /leov1/bites_ma_spain/data/corine_landcover/u2018.....csv 
	--era5_csv /leov1/bites_ma_spain/data/era5/YYYY-MM-DD.csv
	--output /leov1/bites_ma_spain/data/model_input/YYYY-MM-DD.csv 


/leov1/bites_ma_spain/data/model_input/YYYY-MM-DD.csv


R concatDateRange.R
	--from_date 2020-01-01.csv
	--to_date 2025-01-02.csv
	--output /leov1/bites_ma_spain/data/model_input/last_concat.csv


/leov1/bites_ma_spain/data/model_input/last_concat.csv

-------------------------------------------------

MODEL TRAINING (every 15 days)

R train.R
	--data /leov1/bites_ma_spain/data/model_input/last.csv
	--output_dir /leov1/bites_ma_spain/output/brms/train.rds


MODEL PREDICT (1 days)
R predict.R
	--model_rds /leov1/bites_ma_spain/output/brms/train.rds
	--landcover_csv /leov1/bites_ma_spain/data/corine_landcover/u2018.....csv 
	--era_csv /leov1/bites_ma_spain/data/era5/YYYY-MM-DD.csv
	--output /leov1/bites_ma_spain/output/brms/predict/YYYY-MM-DD.csv
