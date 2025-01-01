

## Reproducing results

### Numerical analysis

#### Julia environment

All code was run on Julia v1.10.0 (or v1.9.3 as was the latest available for HPC running of `fig_3_period_over_grid.jl`). A Julia project has been included (`Package.toml`) describing the necessary packages for reproduction. Run the following code to instantiate the necessary packages:

```
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

All necessary Julia code is provided in the `src` subdirectory.

#### Basic figures

The files `fig_2_basic.jl`, `fig_3_bifurcation.jl`, and `fig_5_examples.jl` may now be run. These will produce JLD2 files (HDF5) that we will read from R later to visualise the results.

#### Dynamics over decay rate and seasonality

The figure describing dynamics over varying antibody decay rate and strength of seasonality was produced on the University of Melbourne's SPARTAN high performance computing infrastructure. As such, some code may be specific to the particulars of SPARTAN. The Julia file `fig_3_period_over_grid.jl` performs all analysis, `fig_3_period_over_grid.slurm` defines the SPARTAN job, and `fig_3_period_ver_grid_summ.jl` combines the 8350 outputs into a singular JLD2 file for reading in R.

#### R Environment

All visualisation of the numerical results was performed in R (version 4.2.0), with code provided in the `R` subdirectory. The 