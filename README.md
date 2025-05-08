

## Reproducing results

### Numerical analysis

#### Julia environment

All code was run on Julia v1.11. A Julia project has been included
(`Package.toml`) describing the necessary packages for reproduction. Run the
following code to instantiate the necessary packages:

```
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

All necessary Julia code is provided in the `src` subdirectory.

#### Figures

The files `fig_2_basic.jl`, `fig_3_bifurcation.jl`, `fig_5_examples.jl`, `fig_7_examples.jl` may now be run. These will produce JLD2 files (HDF5) that we will read from R later 
to visualise the results.

The data for other figures can similarly be produced using the `fig_` or `fig_supp_` `.jl` files.

#### Dynamics over decay rate and seasonality

The figure describing dynamics over varying antibody decay rate and strength of 
seasonality was produced on the University of Melbourne's SPARTAN high performance
computing infrastructure. As such, some code may be specific to the particulars 
of SPARTAN. The Julia file `fig_3_period_over_grid.jl` performs all analysis, 
`fig_3_period_over_grid.slurm` defines the SPARTAN job, and `fig_3_period_ver_grid_summ.jl`
combines the outputs into a singular JLD2 (HDF5) file for reading in R.

#### R Environment

All visualisation of the numerical results was performed in R (version 4.2.0),
with code provided in the `R` subdirectory.

Each figures in the main manuscript is produced in the `plot_1_...` files. The
appendix figures are produced in the `plot_supp_...` files.




