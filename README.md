# Code for *Bayesian Unit-level Modeling of Categorical Survey Data with a Longitudinal Design*

This repository contains code associated with **"Bayesian Unit-level Modeling of Categorical Survey Data with a Longitudinal Design"** by **Vedensky, Parker, and Holan**.

## Repository Structure

### `blur/`
Contains all of the model code in **R library format**.

### `code/`
Contains code specific to the manuscript.

- `helper_functions.R`
  Contains helper functions used by the other files.

- `cross_sectional_sim_blur.R`
  Runs the cross-sectional models evaluated in the empirical simulation.

- `longitudinal_sim_blur.R`
  Runs the longitudinal models evaluated in the empirical simulation.

- `process_sim_output.R`
  Generates the figures and tables included in Section 5.2 of the manuscript.

- `data_analysis_blur.R`
  Fits the model to the full HPS data set and produces the figures included in Section 5.3 of the manuscript.

### `data/`
Contains code for processing data, as well as the resulting processed data files.

- `make_basis_functions.R`
  Generates `scaled_basis_functions.RData` and `unscaled_basis_functions.RData`.

- `generate_GAD2_pop.R`
  Generates the dataset used in the data analysis and as the "empirical population" that is subsampled in the simulation study.

- `generate_samples.R`
  Creates the 100 subsamples of the population data used in the empirical simulation. These are stored as `empirical_samples_GAD2.RData`.

- `process_census_table.R`
  Produces `census_tables.csv` for use in `code/data_analysis_blur.R`.


### Note on dataset
The processed data are available as `HPS_empirical_pop_df_GAD2.RData'.
The raw data are too large to be included directly but can be downloaded from `census.gov` with the set of bash commands below then processed with `generate_GAD2_pop.R`.

```
 for i in {1..12}
 do if [ $i -lt 10 ]
 then wget https://www2.census.gov/programs-surveys/demo/datasets/hhp/2020/wk$i/HPS_Week0$i'_PUF_CSV.zip'
 else wget https://www2.census.gov/programs-surveys/demo/datasets/hhp/2020/wk$i/HPS_Week$i'_PUF_CSV.zip'
 fi
 done
 unzip HPS_Week*PUF_CSV.zip
```
