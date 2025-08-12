# Bite prediction - Pipeline HPC

This repository contains the code for running a high-performance computing (HPC) pipeline that predicts the probability of mosquito bites using the [BRSM model](https://github.com/Mosquito-Alert/Bites_MA_Spain).

The pipeline is orchestrated using [Snakemake](https://snakemake.github.io) which automates tasks related to data preparation, model training, and prediction.

## Data sources
The model needs data from the following data sources:
- Mosquito Bite Reports from Mosquito Alert (variable used to train and to be predicted)
- [Sampling effort data](https://github.com/Mosquito-Alert/sampling_effort_data) - from Mosquito Alert.
- [ERA5-Land](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land?tab=overview) - historical climate reanalysis data.
- [CORINE Landcover](https://land.copernicus.eu/en/products/corine-land-cover) – land use/cover information for Europe.

## Data preparation

### Download CORINE Land Cover

Please download the CORINE Land Cover dataset manually and place it in the following path:
`data/corine_landcover/U2018_CLC2018_V2020_20u1.tif`

### Provide vector file (for zonal statistics)

Provide a GeoPackage file containing polygons, with:
- Only one layer
- Two columns: code (unique ID) and name

*IMPORTANT*: Update `config/config.yaml` to reflect the path to this file.

### CDS API Key

To access ERA5-Land data, you'll need an API key from the Copernicus Climate Data Store. Follow instruction at: https://cds.climate.copernicus.eu/how-to-api

## Deployment

### 1.Install Snakemake

We recommend using [Mamba](https://mamba.readthedocs.io/en/latest/) (a faster drop-in replacement for Conda). If you don't have Conda or Mamba installed, consider installing [Miniforge](https://github.com/conda-forge/miniforge).

Install Snakemake, Snakedeploy, and necessary plugins:

```bash
mamba create -c conda-forge -c bioconda --name snakemake snakemake=9.8.1 snakedeploy=0.11.0
```

If you're running on an HPC with SLURM, install additional plugins:
```bash
mamba install -n snakemake -c bioconda snakemake-executor-plugin-slurm=1.5.0 snakemake-storage-plugin-fs=1.1.2
```

Activate the environment:

```bash
conda activate snakemake
```

### 2. Deploy the workflow

Create and move into a project directory:
```bash
mkdir -p path/to/project-workdir
cd path/to/project-workdir
```

Deploy the workflow using Snakedeploy:
```bash
snakedeploy deploy-workflow <URL_TO_THIS_REPO> . --tag <DESIRED_TAG>
```

This will create two directories:
- `workflow/`: contains the deployed Snakemake module
- `config/`: contains configuration files

### 3. Configure workflow

Edit `config/config.yaml` to specify your settings (paths, parameters, etc.) according to your data and environment

### 4. Run workflow

#### Local execution with conda

```bash
snakemake --cores all --sdm conda
```

#### HPC execution with SLURM
Use the provided SLURM profile:
```bash
snakemake --cores all --sdm conda --profile slurm
```

For advanced features such as cluster execution, cloud deployments, and workflow customization, see the [Snakemake documentation](https://snakemake.readthedocs.io/).
