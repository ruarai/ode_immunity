include("dependencies.jl")

using JLD2

k = 32
C = 8.0
R = 1.5
gamma = 0.25
beta = R * gamma
rho = 0.003
b = 2^3
h_0 = 8
c_jump_dist = Normal(2^6, 2^5)

model_params_0 = make_model_parameters(
    k = k, beta = beta, gamma = gamma, C = C, rho = rho,
    b = b, h = h_0, c_jump_dist = c_jump_dist; boosting = "none"
)

plot(model_params_0.p_acq)

ode_sparsity = ode_get_sparsity(model_params_0)

n_inf_0 = 0.01
n_days = 32000



x_h = 2.0 .^ (1.5:0.1:5)

plot(x_h)
plot([model_params_0.c_levels .^ h ./ (model_params_0.b .^ h .+ model_params_0.c_levels .^ h) for h in x_h])



y_fixed_I = zeros(length(x_h))

y_I_sol = zeros(length(x_h), n_days)

Threads.@threads for i in eachindex(x_h)

    model_params = make_model_parameters(
        k = k, beta = beta, gamma = gamma, C = C, rho = rho,
        b = b, h = x_h[i], c_jump_dist = c_jump_dist; boosting = "none"
    )

    # Calculate the fixed point/steady state solution
    # sus_vec, inf_vec = get_steady_state(model_params)
    # y_fixed_I[i] = sum(inf_vec)


    # Calculate a (not-necessarily-stable) solution of the ODE over time
    ode_solution = ode_solve(model_params, n_days, n_inf_0, ode_sparsity)

    for d in 1:n_days
        y_I_sol[i, d] = sum(ode_solution(d)[ode_ix(c_inf, 1:model_params.S, model_params.S)])
    end

    println(i)
end


plot(y_I_sol[8, 30000:end])

heatmap(y_I_sol[:,1:8000])
heatmap(y_I_sol[:,20000:end])

plot(y_I_sol[:,end], legend = false)


