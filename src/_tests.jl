
include("dependencies.jl")

using Test


k = 16
R = 1.5
gamma = 0.25
beta = R * gamma
lambda = 0.002
b = 0.25
m = 40
c_jump_dist = Normal(0.8, 0.05)

model_params = make_model_parameters(
    k = k, beta = beta, gamma = gamma, lambda = lambda,
    b = b, m = m, c_jump_dist = c_jump_dist; boosting = false
)

@testset "Parameters" begin
    # Boosting matrix
    @test all(sum(model_params.M, dims = 1) .== 1)
    @test all(model_params.M .>= 0)

    @test model_params.M[15] == cdf.(c_jump_dist, model_params.c_levels[16]) - cdf.(c_jump_dist, model_params.c_levels[15])

    # Strata variables
    @test all(model_params.c_levels .>= 0) && all(model_params.p_acq .<= 1)
    @test all(model_params.p_acq .>= 0) && all(model_params.p_acq .<= 1)
end


ode_sparsity = ode_get_sparsity(model_params)

n_inf_0 = 0.0001
n_days = 365*20

ode_solution = @time ode_solve(model_params, n_days, n_inf_0, ode_sparsity)

sol_I = zeros(n_days, model_params.S)
sol_S = zeros(n_days, model_params.S)

for d in 1:n_days, i in 1:model_params.S
    sol_S[d, :] = ode_solution(d)[ode_ix(c_sus, 1:model_params.S, model_params.S)]
    sol_I[d, :] = ode_solution(d)[ode_ix(c_inf, 1:model_params.S, model_params.S)]
end


@testset "Solving" begin
    @test all(sol_S .>= 0) && all(sol_S .<= 1)
    @test all(sol_I .>= 0) && all(sol_I .<= 1)

    @test all(sum(sol_S .+ sol_I, dims = 2) .<= 1.0000001)
    @test all(sum(sol_S .+ sol_I, dims = 2) .>= 0.0)
end

model_params_stoch = make_model_parameters(
    k = 16, beta = 0.25 * 1.5, gamma = 0.25, lambda = 0.006,
    b = 0.25, m = 40, c_jump_dist = Normal(0.8, 0.05); boosting = false
)

u_t, inf = ctmc_sim(model_params_stoch, n_days, 10, 0)
sol_S = u_t[ode_ix(c_sus, 1:model_params.S, model_params.S), :]
sol_I = u_t[ode_ix(c_inf, 1:model_params.S, model_params.S), :]

@testset "Stochastic solving" begin
    @test all(sol_S .>= 0)
    @test all(sol_I .>= 0)
    @test all(sum(sol_S .+ sol_I, dims = 1) .== 100000)
end