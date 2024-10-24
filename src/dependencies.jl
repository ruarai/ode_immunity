

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

# Steady state solutions dependencies
#using MultiFloats
using DoubleFloats
using NonlinearSolve
using StaticArrays

# ODE solver dependencies
using DifferentialEquations
using Symbolics
using Plots

# Stochastics/agent-based dependencies
using Random

include("model_globals.jl")
include("steady_state_functions.jl")
include("model_parameters.jl")
include("solver.jl")
include("solve_stochastic.jl")
include("agent_based.jl")
include("dde.jl")


function expand_grid(; kws...)
    names, vals = keys(kws), values(kws)
    return DataFrame(NamedTuple{names}(t) for t in Iterators.product(vals...))
end