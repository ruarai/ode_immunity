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





# Figures for confirmation report


const stoch_pop_size = 30000 # need to re evaluate functions after this
n_days = 1000
model_params = make_model_parameters(
    k = 16, beta = 0.25 * 1.5, gamma = 0.25, lambda = 0.007,
    b = 0.25, m = 40, c_jump_dist = Normal(0.8, 0.5); boosting = false
)
u_t, inf = ctmc_sim(model_params, n_days, 10, 0)

plot(u_t', legend = false)

heatmap(log.(u_t[1:(model_params.S * 2),:] .+ 1))
plot(inf, xlim = (0, 365 * 1.5))


function ctmc_sim_scenario(model_params, n_days, n_inf_0, seed)
    u0 = zeros(Int64, model_params.S * n_compartments)
    u0[ode_ix(c_sus, 1, model_params.S)] = stoch_pop_size - n_inf_0
    u0[ode_ix(c_inf, 1, model_params.S)] = n_inf_0
    
    du = zeros(Int64, model_params.S * n_compartments)
    
    u = copy(u0)
    u_t = zeros(length(du), n_days)
    
    inf = zeros(n_days)
    
    Random.seed!(seed)
    
    for t in 1:n_days
    
        u_t[:,t] = u
        inf[t] = sum(u[ode_ix(c_inf, 1:model_params.S, model_params.S)])

        # If in trivial steady state (all in S_0), do not simulate any more
        if u[ode_ix(c_sus, 1, model_params.S)] == stoch_pop_size
            continue
        end

        beta_mult = 1.0

        if t < 250
            beta_mult = 0.8
        end
    
        for i in 1:100
            ctmc_step!(du, u, model_params, 0.01; beta_mult = beta_mult)
            u = u + du
        end
    end

    return u_t, inf
end

u_t, inf_scen = ctmc_sim_scenario(model_params, n_days, 10, 2)
plot!(inf_scen)


function ctmc_sim_scenario_2(model_params, n_days, n_inf_0, seed)
    u0 = zeros(Int64, model_params.S * n_compartments)
    u0[ode_ix(c_sus, 1, model_params.S)] = stoch_pop_size - n_inf_0
    u0[ode_ix(c_inf, 1, model_params.S)] = n_inf_0
    
    du = zeros(Int64, model_params.S * n_compartments)
    
    u = copy(u0)
    u_t = zeros(length(du), n_days)
    
    inf = zeros(n_days)
    
    Random.seed!(seed)
    
    for t in 1:n_days
    
        u_t[:,t] = u
        inf[t] = sum(u[ode_ix(c_inf, 1:model_params.S, model_params.S)])

        # If in trivial steady state (all in S_0), do not simulate any more
        if u[ode_ix(c_sus, 1, model_params.S)] == stoch_pop_size
            continue
        end

        flow_vacc = 1e-10

        if t > 150 && t < 400
            flow_vacc = 4e-8
        end
    
        for i in 1:100
            ctmc_step!(du, u, model_params, 0.01; flow_vacc = flow_vacc)
            u = u + du
        end
    end

    return u_t, inf
end


u_t, inf_scen_2 = ctmc_sim_scenario_2(model_params, n_days, 10, 2)
plot!(inf_scen_2)


jldsave("data/stochastic_inf_scenarios.jld2"; inf, inf_scen, inf_scen_2)

Plots.pdf("results/stochastic_inf_scenarios.pdf")