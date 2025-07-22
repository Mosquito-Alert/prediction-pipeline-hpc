## Fetch data

1. Download ERA5 for Spain (@epou) -> `./data/era5/spain/YYYY-MM-DD.grib`
2. Convert ERA5 grib to CSV format applying a clipping mask for Spain (@CatuCerecedo)
	```
	R data/era5/grib2csv.R
		--input ./data/era5/spain/YYYY-MM-DD.grib
		--boundary_mask ./data/boundaries/gadm41_ESP_simplified.gpkg
		--boundary_column_id code
		--output ./data/era5/spain/YYYY-MM-DD.csv
	```
3. Fetch bite reports (@epou) -> `./data/bite_reports/YYYY-MM-DD.csv`
4. (Once) Download CORINE landcover -> `./data/corine_landcover/u2018.....tif`
5. (Once) Apply Spain mask to CORINE landcover (@CatuCerecedo) -> `./data/corine_landcover/u2018.....csv`

## Models

### Bite BRSM model0

#### Data preparation
1. Apply lags to ERA5 data (@CatuCerecedo)
	```
	R models/bites_brms0/preprocess/process_era5.R
		--era5_dir ./data/era5/spain/
		--format 'YYYY-MM-DD.csv'
		--output ./outputs/bites_brsm0/processed_era5/YYYY-MM-DD.csv
	```
2. Generate features (@CatuCerecedo)
	```
	R models/bites_brms0/preprocess/generate_features.R
		--bites_csv ./data/bite_reports/YYYY-MM-DD.csv
		--boundary_mask ./data/boundaries/gadm41_ESP_simplified.gpkg 
		--boundary_column_id code 
		--landcover_csv ./data/corine_landcover/u2018.....csv 
		--era5_csv ./outputs/bites_brsm0/processed_era5/YYYY-MM-DD.csv
		--output ./outputs/bites_brsm0/features/YYYY-MM-DD.csv 
	```

#### Train
1. Concat multiple daily features (@CatuCerecedo)
```
	R models/bites_brms0/preprocess/concat_features.R
		--features_dir ./outputs/bites_brsm0/features/
		--format 'YYYY-MM-DD.csv'
		--from_date 2020-01-01
		--to_date 2025-01-02
		--output ./outputs/bites_brsm0/features/last_concat.csv 
```
2. Train (@CatuCerecedo)
```
	R models/bites_brms0/train.R
		--data ./outputs/bites_brsm0/features/last_concat.csv 
		--output_dir ./output/bites_brsm0/model/train.rds
```

#### Predict
1. Predict (@CatuCerecedo)
```
	R predict.R
		--model_rds ./output/bites_brsm0/model/train.rds
		--landcover_csv ./data/corine_landcover/u2018.....csv 
		--era_csv ./data/era5/spain/YYYY-MM-DD.csv
		--output ./output/bites_brsm0/predictions/YYYY-MM-DD.csv
```
