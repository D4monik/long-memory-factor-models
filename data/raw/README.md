# Raw data

The original data used in the paper are from the Global Historical Climatology 
Network monthly temperature dataset, version 4 (GHCN-M v4), published by the 
NOAA National Centers for Environmental Information.

The compressed archive `ghcnm.tavg.v4.0.1.20260330.qcf.zip` contains:

- `ghcnm.tavg.v4.0.1.20260330.qcf.dat`: monthly mean temperatures, quality-controlled and adjusted
- `ghcnm.tavg.v4.0.1.20260330.qcf.inv`: station inventory

These files correspond to version `v4.0.1.20260330`, downloaded on 31 March 2026. 
The archive is included to preserve the exact data version used in the paper, as 
NOAA updates the dataset regularly.

To rebuild the processed panel, delete `data/processed/temperature_panel.csv` and run
`run_all.R`. The archive will be extracted automatically. Alternatively, extract both 
files into `data/raw/` and run `scripts/02_prepare_data.R` directly.

The original dataset is available from [NOAA NCEI](https://www.ncei.noaa.gov/pub/data/ghcn/v4/).
The GHCN-M v4 dataset is listed under [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/) 
in the Data.gov catalog. Please cite the original NOAA dataset when using these data.