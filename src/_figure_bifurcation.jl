
include("src/dependencies.jl")

c_max = 8
k = 4
beta = 1.5
gamma = 0.5
b = 3
m = 4
c_jump_dist = Normal(7, 0.5)

# Temporary model_parameters so sparsity can be generated
model_params_0 = make_model_parameters(
    c_max = c_max, k = k, beta = beta, gamma = gamma, lambda = 0.1,
    b = b, m = m, c_jump_dist = c_jump_dist
)
ode_sparsity = ode_get_sparsity(model_params)

n_inf_0 = 0.01
n_days = 16000



x_lambda = collect(0.01:0.0025:0.5)
#x_lambda = [0.02]

y_fixed_I = zeros(length(x_lambda))

y_I_sol = zeros(length(x_lambda), n_days)

Threads.@threads for i in eachindex(x_lambda)
    println(i)


    model_params = make_model_parameters(
        c_max = c_max, k = k, beta = beta, gamma = gamma, lambda = x_lambda[i],
        b = b, m = m, c_jump_dist = c_jump_dist
    )

    # Calculate the fixed point/steady state solution
    sus_vec, inf_vec = get_steady_state(model_params)
    y_fixed_I[i] = sum(inf_vec)


    # Calculate a (non-stable) solution of the ODE over time
    ode_solution = @time ode_solve(model_params, n_days, n_inf_0, ode_sparsity)

    for d in 1:n_days
        y_I_sol[i, d] = sum(ode_solution(d)[ode_ix(c_inf, 1:model_params.N, model_params.N)])
    end

end



jldsave("data/bifurcations.jld2"; x_lambda, y_fixed_I, y_I_sol)

plot(x_lambda, y_fixed_I)
plot!(x_lambda, maximum(y_I_sol[:, 8000:12000], dims = 2))

plot(y_I_sol[1,:], legend = false)


x_lambda[10]
plot(y_I_sol[10,1:500])


x_lambda[40]
plot(y_I_sol[40,1:500])