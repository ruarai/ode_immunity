

include("dependencies.jl")


using NonlinearSolve, StaticArrays

params = make_model_parameters(
    c_max = 8, k = 4,

    beta = 1.0, sigma = 1.0, gamma = 0.5, lambda = 0.05,

    b = 3, m = 4,

    c_jump_dist = Normal(4, 0.5)
)

B = build_waning_matrix(params)
M = build_immunity_matrix(params)


omega_inv = 1 ./ (1 .- params.p_acq)


u0 = @SVector BigFloat[rand() for x in 1:params.N]

fn_solve(inf_vec, p) = steady_state_and_valid(
    inf_vec, 
    B, M, 
    omega_inv, 
    params.k, params.lambda, params.beta, params.gamma, params.sigma
)

probN = NonlinearProblem(fn_solve, u0)
solver = solve(probN, NewtonRaphson(), abstol = 1e-16)



fn_solve(solver.u, 1)

plot(convert.(Float64, solver.u))