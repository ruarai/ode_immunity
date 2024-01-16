
include("dependencies.jl")



model_params = make_model_parameters(
    c_max = 8, k = 4,

    beta = 1.5, sigma = 1.0, gamma = 0.5, lambda = 0.05,

    b = 4, m = 4,

    c_jump_dist = Normal(6, 0.5)
)


ode_sparsity = ode_get_sparsity(model_params)

n_inf_0 = 0.01
n_days = 3000

ode_solution = @time ode_solve(model_params, n_days, n_inf_0, ode_sparsity)

infected = zeros(300, 32)

for k in 1:model_params.N
    for d in 1:n_days
        infected[d, k] = sum(ode_solution(d)[ode_ix(c_inf, k, model_params.N)])
    end
end


plot(log.(infected))


