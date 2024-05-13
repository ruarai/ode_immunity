include("dependencies.jl")

using JLD2


# Continued periodic behaviour
n_days = 8000
model_params = make_model_parameters(
    k = 16, beta = 0.25 * 1.5, gamma = 0.25, lambda = 0.006,
    b = 0.25, m = 40, c_jump_dist = Normal(0.8, 0.05); boosting = false
)
u_t, inf = ctmc_sim(model_params, n_days, 10, 0)

plot(u_t', legend = false)

heatmap(log.(u_t[1:(model_params.S * 2),:] .+ 1))
plot(inf)

sol_I = zeros(n_days, model_params.S)
sol_S = zeros(n_days, model_params.S)

for d in 1:n_days, i in 1:model_params.S
    sol_S[d, :] = u_t[ode_ix(c_sus, 1:model_params.S, model_params.S), d]
    sol_I[d, :] = u_t[ode_ix(c_inf, 1:model_params.S, model_params.S), d]
end

jldsave("data/paper/stochastic_basic.jld2"; sol_I, sol_S)




# Periodic behaviour around the limit cycle
n_days = 8000
model_params = make_model_parameters(
    k = 16, beta = 0.25 * 1.5, gamma = 0.25, lambda = 0.008,
    b = 0.25, m = 40, c_jump_dist = Normal(0.8, 0.05); boosting = false
)
u_t, inf = ctmc_sim(model_params, n_days, 10, 0)

plot(inf)


sol_I = zeros(n_days, model_params.S)
sol_S = zeros(n_days, model_params.S)

for d in 1:n_days, i in 1:model_params.S
    sol_S[d, :] = u_t[ode_ix(c_sus, 1:model_params.S, model_params.S), d]
    sol_I[d, :] = u_t[ode_ix(c_inf, 1:model_params.S, model_params.S), d]
end

jldsave("data/paper/stochastic_around_equilbrium.jld2"; sol_I, sol_S)


# Extinction

model_params = make_model_parameters(
    k = 16, beta = 0.25 * 1.5, gamma = 0.25, lambda = 0.005,
    b = 0.25, m = 40, c_jump_dist = Normal(0.8, 0.05); boosting = false
)

n_days = 4000
n_sims = 400

inf_sims = zeros(n_days, n_sims)

Threads.@threads for i in 1:n_sims
    println(i)
    u_t, inf = ctmc_sim(model_params, n_days, 10, i)
    inf_sims[:, i] = inf
end


jldsave("data/paper/stochastic_extinction_sims.jld2"; inf_sims)

heatmap(log.(inf_sims .+ 1)')
plot(inf_sims, legend = false, lc = "black")



u_t, inf = ctmc_sim(model_params, n_days, 10, 0)
plot(inf)


sol_I = zeros(n_days, model_params.S)
sol_S = zeros(n_days, model_params.S)

for d in 1:n_days, i in 1:model_params.S
    sol_S[d, :] = u_t[ode_ix(c_sus, 1:model_params.S, model_params.S), d]
    sol_I[d, :] = u_t[ode_ix(c_inf, 1:model_params.S, model_params.S), d]
end

jldsave("data/paper/stochastic_extinction.jld2"; sol_I, sol_S)


