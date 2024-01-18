

include("src/dependencies.jl")



model_params = make_model_parameters(
    c_max = 8, k = 4,

    beta = 1.0, sigma = 1.0, gamma = 0.5, lambda = 0.05,

    b = 4, m = 4,

    c_jump_dist = Normal(6, 0.5)
)

solver = get_steady_state(model_params)

fn_solve(solver.u, 1)