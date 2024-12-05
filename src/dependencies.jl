
using Pkg

if dirname(Base.active_project()) != pwd()
    Pkg.activate(".")
end


if VERSION != v"1.10.0"
    println("Julia version has changed. Is this correct?")
end


# Base dependencies
using Distributions
using LinearAlgebra
using ProgressMeter
using NaNMath
using JLD2
using Dates

# Steady state solutions dependencies
#using MultiFloats
using DoubleFloats
using NonlinearSolve
using StaticArrays

# ODE solver dependencies
using DifferentialEquations
using Symbolics
using Plots
using ChaosTools

# Stochastics/agent-based dependencies
using Random

include("model_globals.jl")
include("steady_state_functions.jl")
include("model_parameters.jl")
include("solver.jl")

include("model_no_boosting.jl")

include("helpers.jl")


function get_jobs(arg_ix, n_array, n_jobs)
    N = n_jobs
    M = n_array
    n = ceil(Int, convert(Float64, N) / M)

    return [((i - 1) * n + 1):(min(i * n, N)) for i in 1:M][arg_ix]
end

function expand_grid(; kws...)
    names, vals = keys(kws), values(kws)
    return DataFrame(NamedTuple{names}(t) for t in Iterators.product(vals...))
end


println("$(now()) -- dependencies loaded")