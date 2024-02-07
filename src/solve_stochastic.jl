
const stoch_pop_size = 100000

function ctmc_step!(du, u, model_params, dt)
    k = model_params.k

    for i in 1:(k * n_compartments)
        du[i] = 0
    end

    counts_decay = @views rand.(Binomial.(u[ode_ix(1, 2:k, k)], model_params.lambda .* k .* dt))

    du[ode_ix(1, 2:k, k)] = -counts_decay
    du[ode_ix(1, 1:(k - 1), k)] .+= counts_decay

    flow_sus_to_inf = @views sum(u[ode_ix(c_inf, 1:k, k)]) .*
        model_params.beta .* 
        (1 .- model_params.p_acq) ./
        stoch_pop_size

    count_sus_to_inf = @views rand.(Binomial.(u[ode_ix(c_sus, 1:k, k)], -expm1.(-flow_sus_to_inf .* dt)))
        

    du[ode_ix(c_sus, 1:k, k)] .+= -count_sus_to_inf
    du[ode_ix(c_inf, 1:k, k)] .+= count_sus_to_inf 

    du[ode_ix(c_count, 1:k, k)] .+= count_sus_to_inf 

    for j in 1:k
        count_inf_depart = @views rand(Binomial(u[ode_ix(c_inf, j, k)], -expm1.(-model_params.gamma .* dt)))
        
        du[ode_ix(c_inf, j, k)] -= count_inf_depart

        count_sus_arrive = @views rand(Multinomial(count_inf_depart, model_params.M[:,j], check_args = false))

        for i in 1:k
            du[ode_ix(c_sus, i, k)] += count_sus_arrive[i]
        end

    end
end

function ctmc_sim(model_params, n_days, n_inf_0, seed)
    u0 = zeros(Int64, model_params.k * n_compartments)
    u0[ode_ix(c_sus, 1, model_params.k)] = stoch_pop_size - n_inf_0
    u0[ode_ix(c_inf, 1, model_params.k)] = n_inf_0
    
    du = zeros(Int64, model_params.k * n_compartments)
    
    u = copy(u0)
    u_t = zeros(length(du), n_days)
    
    inf = zeros(n_days)
    
    Random.seed!(seed)
    
    for t in 1:n_days
    
        u_t[:,t] = u
        inf[t] = sum(u[ode_ix(c_inf, 1:model_params.k, model_params.k)])

        # If in trivial steady state (all in S_0), do not simulate any more
        if u[ode_ix(c_sus, 1, model_params.k)] == stoch_pop_size
            continue
        end
    
        for i in 1:100
            ctmc_step!(du, u, model_params, 0.01)
            u = u + du
        end
    end

    return u_t, inf
end
