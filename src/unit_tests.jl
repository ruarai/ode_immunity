
include("dependencies.jl")
using Test


model_params = make_model_parameters(
    k = baseline_k, beta = baseline_beta, gamma = baseline_gamma,
    a = baseline_a, r = baseline_r,
    b = baseline_b, h = baseline_h, c_jump_dist = baseline_c_jump_dist
)

@testset "Parameters" begin
    # Boosting matrix
    @test all(sum(model_params.M, dims = 1) .>= 1 - 1e-10)
    @test all(model_params.M .>= 0)
    
    @test all(sum(model_params.p_trans) .>= 1 - 1e-10)
    @test all(sum(model_params.p_trans) .<= 1 + 1e-10)

    # Strata variables
    @test all(model_params.c_levels .>= 0) && all(model_params.p_acq .<= 1)
    @test all(model_params.p_acq .>= 0) && all(model_params.p_acq .<= 1)
end


ode_solution = @time ode_solve(model_params, n_days, n_inf_0)

t = 1:n_days
sus, inf, inc = get_results(ode_solution, t, model_params)

@testset "Model result sanity" begin
    @test all(sus .>= 0) && all(sus .<= 1)
    @test all(inf .>= 0) && all(inf .<= 1)
    @test all(inc .>= 0) && all(inc .<= 1)

    @test all(sum(sus, dims = 2) .+ inf .<= 1.0000001)
end