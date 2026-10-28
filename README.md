# Supplement B: Code and Data

This supplement provides a local version of the Bagel software used in the paper, together with the code that can reproduce the motivating example, simulation studies, and Azure LLM inference traces analysis. We additionally provide parameter settings, calibrated thresholds, and numerical results, so the results can be reproduced directly from saved outputs to avoid heavy computation. Parralleliztation is implemented in every script, so changing  `N_cores` to your hardware can reduce the running time and avoid crashes.

## Installation

The local package require R and a C++20-compatible compiler. On Windows, Rtools is required to compile the source code. In the paper, we use the local package in `bagel/bagelR` for simulation. To reproduce the reported analyses, install the package from root.
```r
install.packages(c("Rcpp", "RcppEigen", "BH", "dplyr", "ggplot2", "ggrepel", "paletteer"))
install.packages("bagel/bagelR", repos = NULL, type = "source")
```
The updated Bagel package can be installed by
```r
if(system.file(package = "bagelR") != "")
    {
      remove.packages("bagelR")
    }
library("remotes")
install_github("grosed/bagel/bagelR",force=TRUE)
```
To make sure obtain consistent results as the paper, please use the local version as it is frozen here.

## Motivating example
The resulting figures are shown in the Figure 2 and Figure 3 in main paper, they can be generated as:
```r
source("simulation/motivation_example.R")
```

## Generating the simulation results

### Reproducing the table about accuracy of the posterior approximation
The results are shown in the Table 1 in main paper. We consider the following scenarios.
| Scenario | Script in `simulation/` |
| --- | --- |
| Mean change with known variance | `example1_varknown_calibration.R` |
| Mean change with unknown variance | `example1_varunknown_calibration.R` |
| Continuous slope change | `example2_calibration_conti.R` |
| Discontinuous slope change | `example2_calibration_dis.R` |
| Seasonal change | `example_seasonality_calibration.R` |

Take mean change with known variance as an example, we can reproduce the result by executing:

```r
CALIB_NO_RUN <- TRUE
source("simulation/example1_varknown_calibration.R")
SKIP_EXISTING <- FALSE
run_calibration()
```

`CALIB_NO_RUN <- TRUE` loads the function definitions without starting the experiment immediately. Setting `SKIP_EXISTING <- FALSE` rerun the simulation and overwrite exisiting results. `REUSE_THRESH <- TRUE` is set within the script. Changing it to false will restart the threshold calibration, which is computationally heavy. The command line is similiar for other scenarios. Results for each scenario will be written to `simulation/results_calibration/`.

After all five scenarios have completed, Table 1 can be generated with:

```r
source("simulation/results_calibration/calibration_summary.R")
```

### Reproducing the table about accuracy under online setting
The results are shown in the Table 3 in main paper and Table 1 in supplement A. Here we consider the following five scenarios. Each script contains two parts: one is to generate the simulation results for exact approach, Bagel, LP-P and LP-M, and another part is to generate the simulation results of exact approach that stops at the stopping time of each approximated approach for each replicate. 
| Scenario | Script in `simulation/` | 
| --- | --- | 
| Mean change with known variance | `example1_varknown.R` | 
| Mean change with unknown variance | `example1_varunknown.R` | 
| Continuous slope change | `example2_conti.R` | 
| Discontinuous slope change | `example2_dis.R` | 
| Seasonal change | `example_seasonality.R` | 

Take mean change with known variance as an example. Within the script, uncommenting `run_example1_varknown()`and execute:
```r
source("simulation/example1_varknown.R")
```
This will run the part 1 and part 2 together and results are written into the folder named `results_example1_varknown`. Similiarly, `REUSE_THRESH <- TRUE` can be set within the script to reuse the provided threhsold. One could also skip the first part by commenting `run_example1_varknown()` and execute the same command line.

After all five scenarios have completed, Table 3 in the main paper and Table 1 in supplement A can be generated with: 
```r
source("simulation/example2_summary.R")
```
Note: this isimulation is comvcomputationally heavy, we suggest use the existing threshold to reproduce the results.

## Empirical runtime

The resulting figures are shown in the Figure 4 in main paper and Figure 2 in supplement A, they can be generated separately with:

```r
source("simulation/speed.R")
```

The generated running speed will not be exactly the same as the paper shows as the runtime depends on the hardware of computers. However, the main trend should be the same.

## Azure LLM inference trace data and analysis

The results are shown in the real data example in main paper. Azure LLM trace data is downloaded from  Jovan Stojkovic, Chaojie Zhang, Íñigo Goiri, Josep Torrellas, Esha Choukse. ["DynamoLLM: Designing LLM Inference Clusters for Performance and Energy Efficiency", in Proceedings of the International Symposium on High-Performance Computer Architecture (HPCA 2025). IEEE, Las Vegas, USA, 2025](https://ieeexplore.ieee.org/abstract/document/10946802/).

For more details and access to the data, see https://github.com/Azure/AzurePublicDataset/blob/master/AzureLLMInferenceDataset2024.md

The results are shown in the real data example section in main paper. Script `training_loss/bin_tokens.py` aggregates the Azure inference trace into five-minute intervals. It requires Python, pandas, to fetch the source data from the Github repo. The processed observations are saved and provided in `training_loss/InferenceTrace_code_5min.csv`.

The parameter settings used in the paper are also provided in `training_loss/code_gen_slope_paper_parameters.txt`. To use these saved parameters instead of tuning threshold from beginning, set `use_saved_parameters <- TRUE` within the script. 

To reproduce the results in paper, execute:
```r
source("training_loss/code_gen_slope_dynamic50.R")
```
