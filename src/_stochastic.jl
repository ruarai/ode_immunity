include("dependencies.jl")

using Random
using JLD2

const stoch_pop_size = 100000

function ctmc_step!(du, u, model_params, dt)
    N = model_params.N

    for i in 1:(N * n_compartments)
        du[i] = 0
    end

    counts_decay = @views rand.(Binomial.(u[ode_ix(1, 2:N, N)], model_params.lambda .* model_params.k .* dt))

    du[ode_ix(1, 2:N, N)] = -counts_decay
    du[ode_ix(1, 1:(N - 1), N)] .+= counts_decay

    flow_sus_to_inf = @views sum(u[ode_ix(c_inf, 1:N, N)]) .*
        model_params.beta .* 
        (1 .- model_params.p_acq) ./
        stoch_pop_size

    count_sus_to_inf = @views rand.(Binomial.(u[ode_ix(c_sus, 1:N, N)], -expm1.(-flow_sus_to_inf .* dt)))
        

    du[ode_ix(c_sus, 1:N, N)] .+= -count_sus_to_inf
    du[ode_ix(c_inf, 1:N, N)] .+= count_sus_to_inf 

    du[ode_ix(c_count, 1:N, N)] .+= count_sus_to_inf 

    for j in 1:N
        count_inf_depart = @views rand(Binomial(u[ode_ix(c_inf, j, N)], -expm1.(-model_params.gamma .* dt)))
        
        du[ode_ix(c_inf, j, N)] -= count_inf_depart

        count_sus_arrive = @views rand(Multinomial(count_inf_depart, model_params.M[:,j], check_args = false))

        for i in 1:N
            du[ode_ix(c_sus, i, N)] += count_sus_arrive[i]
        end

    end
end

function ctmc_sim(model_params, n_days, n_inf_0, seed)
    u0 = zeros(Int64, model_params.N * n_compartments)
    u0[ode_ix(c_sus, 1, model_params.N)] = stoch_pop_size - n_inf_0
    u0[ode_ix(c_inf, 1, model_params.N)] = n_inf_0
    
    du = zeros(Int64, model_params.N * n_compartments)
    
    u = copy(u0)
    u_t = zeros(length(du), n_days)
    
    inf = zeros(n_days)
    
    Random.seed!(seed)
    
    for t in 1:n_days
    
        u_t[:,t] = u
        inf[t] = sum(u[ode_ix(c_inf, 1:model_params.N, model_params.N)])

        # If in trivial steady state (all in S_0), do not simulate any more
        if u[ode_ix(c_sus, 1, model_params.N)] == stoch_pop_size
            continue
        end
    
        for i in 1:100
            ctmc_step!(du, u, model_params, 0.01)
            u = u + du
        end
    end

    return u_t, inf
end


# Continued periodic behaviour
n_days = 8000
model_params = make_model_parameters(
    k = 16, beta = 0.25 * 1.5, gamma = 0.25, lambda = 0.006,
    b = 0.25, m = 40, c_jump_dist = Normal(0.8, 0.05); boosting = false
)
u_t, inf = ctmc_sim(model_params, n_days, 10, 0)

plot(u_t', legend = false)

heatmap(log.(u_t[1:(model_params.N * 2),:] .+ 1))
plot(inf)

sol_I = zeros(n_days, model_params.N)
sol_S = zeros(n_days, model_params.N)

for d in 1:n_days, i in 1:model_params.N
    sol_S[d, :] = u_t[ode_ix(c_sus, 1:model_params.N, model_params.N), d]
    sol_I[d, :] = u_t[ode_ix(c_inf, 1:model_params.N, model_params.N), d]
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


sol_I = zeros(n_days, model_params.N)
sol_S = zeros(n_days, model_params.N)

for d in 1:n_days, i in 1:model_params.N
    sol_S[d, :] = u_t[ode_ix(c_sus, 1:model_params.N, model_params.N), d]
    sol_I[d, :] = u_t[ode_ix(c_inf, 1:model_params.N, model_params.N), d]
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


sol_I = zeros(n_days, model_params.N)
sol_S = zeros(n_days, model_params.N)

for d in 1:n_days, i in 1:model_params.N
    sol_S[d, :] = u_t[ode_ix(c_sus, 1:model_params.N, model_params.N), d]
    sol_I[d, :] = u_t[ode_ix(c_inf, 1:model_params.N, model_params.N), d]
end

jldsave("data/paper/stochastic_extinction.jld2"; sol_I, sol_S)


