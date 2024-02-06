include("dependencies.jl")

using Random
using JLD2

function ctmc_step!(du, u, model_params, dt)
    N = model_params.N

    for i in 1:(N * n_compartments)
        du[i] = 0
    end

    counts_decay = rand.(Binomial.(u[ode_ix(1, 2:N, N)], model_params.lambda .* model_params.k .* dt))

    du[ode_ix(1, 2:N, N)] = -counts_decay
    du[ode_ix(1, 1:(N - 1), N)] .+= counts_decay

    flow_sus_to_inf = @views sum(u[ode_ix(c_inf, 1:N, N)]) .*
        model_params.beta .* 
        (1 .- model_params.p_acq) ./
        100000

    count_sus_to_inf = rand.(Binomial.(u[ode_ix(c_sus, 1:N, N)], -expm1.(-flow_sus_to_inf .* dt)))
        

    mat_inf_to_sus = -expm1.(-dt .* model_params.gamma .* model_params.M)

    # Seems like it would be slow
    cols_inf_to_sus = [vcat(mat_inf_to_sus[:,j], 1 - sum(mat_inf_to_sus[:,j])) for j in axes(mat_inf_to_sus, 2)]

    dist_inf_to_sus = Multinomial.(u[ode_ix(c_inf, 1:N, N)], cols_inf_to_sus)
    count_inf_to_sus = rand.(dist_inf_to_sus)
    count_inf_to_sus = mapreduce(permutedims, vcat, count_inf_to_sus)[:,1:(end - 1)]

    du[ode_ix(c_sus, 1:N, N)] .+= sum(count_inf_to_sus, dims = 1)[1,:] - count_sus_to_inf
    du[ode_ix(c_inf, 1:N, N)] .+= count_sus_to_inf - sum(count_inf_to_sus, dims = 2)[:,1]

    du[ode_ix(c_count, 1:N, N)] .+= count_sus_to_inf 
end

function ctmc_sim(model_params, n_days, n_inf_0, seed)


    u0 = zeros(Int64, model_params.N * n_compartments)
    u0[ode_ix(c_sus, 1, model_params.N)] = 100000 - n_inf_0
    u0[ode_ix(c_inf, 1, model_params.N)] = n_inf_0
    
    du = zeros(Int64, model_params.N * n_compartments)
    
    u = copy(u0)
    u_t = zeros(length(du), n_days)
    
    inf = zeros(n_days)
    
    Random.seed!(seed)
    
    for t in 1:n_days
    
        u_t[:,t] = u
        inf[t] = sum(u[ode_ix(c_inf, 1:model_params.N, model_params.N)])

        if inf[t] == 0
            break
        end
    
        for i in 1:100
            ctmc_step!(du, u, model_params, 0.01)
            u = u + du
        end
    end

    return u_t, inf
end


# Continued periodic behaviour
# n_days = 20000
# model_params = make_model_parameters(
#     c_max = 8, k = 2, beta = 0.25 * 1.5, gamma = 0.25, lambda = 0.03,
#     b = 2, m = 5, c_jump_dist = Normal(4.5, 0.5)
# )
# u_t, inf = ctmc_sim(model_params, n_days, 10, 0)

# heatmap(log.(u_t[1:(model_params.N * 2),:] .+ 1))
# plot(inf)

# sol_I = zeros(n_days, model_params.N)
# sol_S = zeros(n_days, model_params.N)

# for d in 1:n_days, i in 1:model_params.N
#     sol_S[d, :] = u_t[ode_ix(c_sus, 1:model_params.N, model_params.N), d]
#     sol_I[d, :] = u_t[ode_ix(c_inf, 1:model_params.N, model_params.N), d]
# end

# jldsave("data/paper/stochastic_basic.jld2"; sol_I, sol_S)




# Periodic behaviour around the limit cycle
# n_days = 20000
# model_params = make_model_parameters(
#     c_max = 8, k = 2, beta = 0.25 * 1.5, gamma = 0.25, lambda = 0.04,
#     b = 2, m = 5, c_jump_dist = Normal(4.5, 0.5)
# )
# u_t, inf = ctmc_sim(model_params, n_days, 10, 0)

# plot(inf)


# sol_I = zeros(n_days, model_params.N)
# sol_S = zeros(n_days, model_params.N)

# for d in 1:n_days, i in 1:model_params.N
#     sol_S[d, :] = u_t[ode_ix(c_sus, 1:model_params.N, model_params.N), d]
#     sol_I[d, :] = u_t[ode_ix(c_inf, 1:model_params.N, model_params.N), d]
# end

# jldsave("data/paper/stochastic_around_equilbrium.jld2"; sol_I, sol_S)


# Extinction

model_params = make_model_parameters(
    c_max = 8, k = 2, beta = 0.25 * 1.5, gamma = 0.25, lambda = 0.03,
    b = 2, m = 5, c_jump_dist = Normal(4.5, 0.5)
)

n_days = 5000
n_sims = 500

inf_sims = zeros(n_days, n_sims)

Threads.@threads for i in 1:n_sims
    println(i)
    u_t, inf = ctmc_sim(model_params, n_days, 10, i)
    inf_sims[:, i] = inf

end


jldsave("data/paper/stochastic_extinction_sims.jld2"; inf_sims)

heatmap(log.(inf_sims .+ 1)')
plot(inf_sims, legend = false, lc = "black")



u_t, inf = ctmc_sim(model_params, n_days, 10, 12)
plot(inf)


sol_I = zeros(n_days, model_params.N)
sol_S = zeros(n_days, model_params.N)

for d in 1:n_days, i in 1:model_params.N
    sol_S[d, :] = u_t[ode_ix(c_sus, 1:model_params.N, model_params.N), d]
    sol_I[d, :] = u_t[ode_ix(c_inf, 1:model_params.N, model_params.N), d]
end

jldsave("data/paper/stochastic_extinction.jld2"; sol_I, sol_S)