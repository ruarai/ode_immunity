using Pkg
Pkg.activate(".")

# Base dependencies
using Distributions
using LinearAlgebra

# Steady state solutions dependencies
#using MultiFloats
using DoubleFloats
using NonlinearSolve
using StaticArrays

# ODE solver dependencies
using DifferentialEquations
using Symbolics
using Plots

include("model_globals.jl")
include("steady_state_functions.jl")
include("model_parameters.jl")
include("solver.jl")


function expand_grid(; kws...)
    names, vals = keys(kws), values(kws)
    return DataFrame(NamedTuple{names}(t) for t in Iterators.product(vals...))
end