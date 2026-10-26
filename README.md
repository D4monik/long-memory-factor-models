# Estimation and Model Selection in Factor Models with Long Memory

Replication code for

> Ammon, D., Hartl, T. and Tschernig, R. (2026). Estimation and Model Selection in Factor Models with Long Memory. Working paper.

The repository reproduces the Monte Carlo study (Section 4) and the application 
to European temperatures (Section 5).

## Repository structure

```
├── run_all.R                   runs everything
├── R/                          functions
│   ├── config.R                paths and settings
│   ├── factor_selection.R      Bai-Ng criterion and selection of Alessi et al. (2010)
│   ├── simulation.R            data-generating process and threshold search
│   ├── estimation.R            low-frequency transformation, PCA, warming measures
│   ├── lowfreq_memory.R        low-frequency likelihood for the memory parameter
│   └── plots.R                 figures
├── scripts/
│   ├── 01_simulation.R         Section 4, Table 1
│   ├── 02_prepare_data.R       Section 5.1: balanced station panel
│   ├── 03_levels.R             Sections 5.1-5.2: factor model in levels
│   ├── 04_memory.R             Section 5.3: memory of the factors
│   └── 05_differences.R        Section 5.4: factor model after first differencing
├── data/
│   ├── raw/                    original GHCN-M data (compressed archive)
│   └── processed/              balanced panel and station metadata
└── output/                     tables and figures (created by the scripts)
```

## Requirements

The results were computed with R 4.5.3. The estimation uses only base R. The
following packages are used for the figures and, in `02_prepare_data.R`, for
reading the raw data (versions used for the paper in parentheses):
ggplot2 (4.0.2), sf (1.1-3), rnaturalearth (1.2.0), rnaturalearthdata (1.0.0),
scales (1.4.0), readr (2.2.0), dplyr (1.2.0), stringr (1.6.0), tidyr (1.3.2).

**Quick installation.** Install the current versions with

```r
install.packages(c("ggplot2", "sf", "rnaturalearth", "rnaturalearthdata", "scales",
                   "readr", "dplyr", "stringr", "tidyr"))
```

The project uses renv, which may report that the library is out of sync with
`renv.lock`. This message can be ignored.

**Exact versions.** To install the exact versions recorded in `renv.lock`, open
`long-memory-factor-models.Rproj` in RStudio with R 4.5.x (some recorded versions
cannot be compiled under R 4.6 or later) and run

```r
renv::restore()
```

On Windows, twelve packages (Rcpp, S7, bit64, curl, dplyr, glue, magrittr,
openssl, purrr, s2, stringi, vctrs) are compiled from source, which takes
several minutes and requires
[Rtools45](https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html).
If renv asks whether to install the latest available versions of packages that
failed, answer `n` and resolve the error first.
  
## Data

The temperatures are from the Global Historical Climatology Network monthly dataset, version 4 (GHCN-M v4), mean temperatures (TAVG), quality-controlled and adjusted, file version `v4.0.1.20260330`, downloaded from the NOAA National Centers for Environmental Information on 31 March 2026. NOAA updates the dataset daily and does not retain these versions long-term, so the processed panel used in the paper is included in `data/processed/`:

- `temperature_panel.csv`: monthly mean temperatures in °C, January 1960 to December 2025, one row per month and one column per station (308 stations)
- `station_metadata.csv`: station identifier, country code, name, latitude, longitude and elevation

The original GHCN-M v4 source files used in the paper are included as a compressed archive in `data/raw/`. The processed panel is derived from these files. The GHCN-M v4 dataset is listed under the Creative Commons CC0 1.0 license in the
[Data.gov catalog](https://catalog.data.gov/dataset/global-historical-climatology-network-monthly-temperature-version-4).
Please cite the dataset and its description below when using these data. NOAA and NCEI provide the data without warranty.

> Menne, M. J., Gleason, B. E., Lawrimore, J., Rennie, J. and Williams, C. N. (2017). Global Historical Climatology Network – Monthly Temperature, Version 4 (TAVG, quality-controlled and adjusted, file version v4.0.1.20260330). NOAA National Centers for Environmental Information. doi:10.7289/V5XW4GTH. Accessed 31 March 2026.
>
> Menne, M. J., Williams, C. N., Gleason, B. E., Rennie, J. J. and Lawrimore, J. H. (2018). The Global Historical Climatology Network monthly temperature dataset, version 4. *Journal of Climate* 31(24): 9835–9854.

## Reproducing the results

From the repository root, run

```
Rscript run_all.R
```

or run the scripts in `scripts/` individually, in the order of their numbers.

- **Monte Carlo study.** For each design and panel size, `01_simulation.R` 
searches by bisection for the idiosyncratic variance at which $IC_{p2}$ selects 
the true number of factors in half of the 1,000 replications. 
The replications run in parallel on all but one core, and replication $r$ uses 
the seed $42 + r$, so that the results do not depend on the number of cores. 
The study takes many hours. Finished cells are saved in `output/simulation_thresholds.rds`,
and an interrupted run continues where it stopped. Set `RUN_SIMULATION <- FALSE` 
in `run_all.R` to skip it.

- **Application.** The scripts start from the processed panel. To rebuild the 
processed panel from the original NOAA data, delete `data/processed/temperature_panel.csv`
before running the scripts again. Most of the run time is spent on the factor-number
selection, which uses 15,400 random subpanels in Sections 5.2 and 5.4.

| Result in the paper | Script | Output (in `output/`) |
|---|---|---|
| Table 1 | `01_simulation.R` | `tables/table1_thresholds.tex` |
| Figure 1 | `03_levels.R` | `figures/fig1a_temperatures.pdf`, `figures/fig1b_lowfreq_components.pdf` |
| Figure 2 | `03_levels.R` | `figures/fig2a_abc_levels.pdf`, `figures/fig2b_eigenvalues_levels.pdf` |
| Figure 3 | `03_levels.R` | `figures/fig3a_factor1.pdf`, `figures/fig3b_factor2.pdf` |
| Figure 4 | `03_levels.R` | `figures/fig4a_common_component.pdf`, `figures/fig4b_factor_contributions.pdf` |
| Figure 5 | `03_levels.R` | `figures/fig5a_loadings_factor1.pdf`, `figures/fig5b_loadings_factor2.pdf` |
| Figure 6 | `03_levels.R` | `figures/fig6a_station_warming.pdf`, `figures/fig6b_regional_warming.pdf` |
| Figure 7 | `04_memory.R` | `figures/fig7_memory.pdf`, `tables/memory_estimates.csv` |
| Figure 8 | `05_differences.R` | `figures/fig8a_abc_differences.pdf`, `figures/fig8b_levels_vs_differences.pdf` |
| Table 2 | `05_differences.R` | `tables/table2_regional_warming.tex` |
| Figure 9 | `05_differences.R` | `figures/fig9_difference_map.pdf` |

The numbers reported in the text of Section 5 are printed to the console by `03_levels.R` and `05_differences.R`.

## License

The R code is released under the MIT License (see `LICENSE`).

The original NOAA GHCN-M v4 data are listed under Creative Commons CC0 1.0. 
The derived data in `data/processed/` are also made available under CC0 1.0.

When using the data, please cite the original NOAA dataset and the accompanying 
publication listed in the Data section.