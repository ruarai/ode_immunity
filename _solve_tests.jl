
include("dependencies.jl")



ode_sparsity = ode_get_sparsity(model_params)





n_days = 10000
n_inf_0 = 0.01
x_lambda = 0.0:0.1:0.8
x_beta = 1.0

infected = zeros(length(x_lambda), length(x_beta), n_days)

for i in eachindex(x_lambda), j in eachindex(x_beta)
    println(i)

    model_params = make_model_parameters(
        c_max = 8, k = 4,
    
        beta = x_beta[j], sigma = 1.0, gamma = 0.5, lambda = x_lambda[i],
    
        b = 4, m = 4,
    
        c_jump_dist = Normal(6, 0.5)
    )

    ode_solution = @time ode_solve(model_params, n_days, n_inf_0, ode_sparsity)

    for d in 1:n_days
        infected[i, j, d] = sum(ode_solution(d)[ode_ix(c_inf, 1:model_params.N, model_params.N)])
    end
end


plot(infected[:,1,:]')
