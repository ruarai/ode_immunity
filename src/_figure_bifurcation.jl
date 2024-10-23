include("dependencies.jl")

using JLD2

k = 32
R = 1.5
gamma = 0.25
beta = R * gamma
b = 0.25
m = 40
c_jump_dist = Normal(0.8, 0.05)

# Temporary model_parameters so sparsity can be generated
model_params_0 = make_model_parameters(
    k = k, beta = beta, gamma = gamma, rho = 0.1,
    b = b, m = m, c_jump_dist = c_jump_dist; boosting = false
)
ode_sparsity = ode_get_sparsity(model_params_0)

n_inf_0 = 0.01
n_days = 32000



x_rho = collect(0.0:0.0001:0.015)

y_fixed_I = zeros(length(x_rho))

y_I_sol = zeros(length(x_rho), n_days)

Threads.@threads for i in eachindex(x_rho)
    model_params = make_model_parameters(
        k = k, beta = beta, gamma = gamma, rho = x_rho[i],
        b = b, m = m, c_jump_dist = c_jump_dist, boosting = false
    )

    # Calculate the fixed point/steady state solution
    sus_vec, inf_vec = get_steady_state(model_params)
    y_fixed_I[i] = sum(inf_vec)


    # Calculate a (not-necessarily-stable) solution of the ODE over time
    ode_solution = ode_solve(model_params, n_days, n_inf_0, ode_sparsity)

    for d in 1:n_days
        y_I_sol[i, d] = sum(ode_solution(d)[ode_ix(c_inf, 1:model_params.S, model_params.S)])
    end

    println(i)
end


jldsave("data/paper/bifurcations.jld2"; x_rho, y_fixed_I, y_I_sol)

plot(x_rho, y_fixed_I, xlim = (0, 0.01))
plot!(x_rho, maximum(y_I_sol[:, 30000:end], dims = 2))



heatmap(y_I_sol[:,1:8000])
heatmap(y_I_sol[:,20000:end])

plot(y_I_sol[1,:], legend = false)