
include("dependencies.jl")
using Test


model_params = make_model_parameters(
    k = baseline_k, beta = baseline_beta, gamma = baseline_gamma,
    C = baseline_C, r = baseline_r,
    b = baseline_b, h = baseline_h, c_jump_dist = baseline_c_jump_dist
)

@testset "Parameters" begin
    # Boosting matrix
    @test all(sum(model_params.M, dims = 1) .>= 1-1e-10)
    @test all(model_params.M .>= 0)

    # @test model_params.M[15] == cdf.(c_jump_dist, log10(model_params.c_levels[16])) - cdf.(c_jump_dist, log10(model_params.c_levels[15]))

    # Strata variables
    @test all(model_params.c_levels .>= 0) && all(model_params.p_acq .<= 1)
    @test all(model_params.p_acq .>= 0) && all(model_params.p_acq .<= 1)
end


ode_sparsity = ode_get_sparsity(model_params)

n_inf_0 = 0.0001
ode_solution = @time ode_solve(model_params, n_days, n_inf_0, ode_sparsity)

t = 1:n_days
sol_I[:, :] = ode_solution(t)[ode_ix(c_sus, 1:model_params.S, model_params.S), :]'
sol_S[:, :] = ode_solution(t)[ode_ix(c_inf, 1:model_params.S, model_params.S), :]'

@testset "Solving" begin
    @test all(sol_S .>= 0) && all(sol_S .<= 1)
    @test all(sol_I .>= 0) && all(sol_I .<= 1)

    @test all(sum(sol_S .+ sol_I, dims = 2) .<= 1.0000001)
    @test all(sum(sol_S .+ sol_I, dims = 2) .>= 0.0)
end