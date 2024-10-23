
const stoch_pop_size = 100000

function ctmc_step!(du, u, model_params, dt; beta_mult = 1.0, flow_vacc = 0.0)
    k = model_params.S

    for i in 1:(k * n_compartments)
        du[i] = 0
    end

    counts_decay = @views rand.(Binomial.(u[ode_ix(1, 2:k, k)], model_params.wane_transition_rate .* dt))

    du[ode_ix(1, 2:k, k)] = -counts_decay
    du[ode_ix(1, 1:(k - 1), k)] .+= counts_decay

    flow_sus_to_inf = @views sum(u[ode_ix(c_inf, 1:k, k)]) .*
        model_params.beta .* 
        beta_mult .*
        (1 .- model_params.p_acq) ./
        stoch_pop_size

    count_sus_to_inf = @views rand.(Binomial.(u[ode_ix(c_sus, 1:k, k)], -expm1.(-flow_sus_to_inf .* dt)))

    sus_not_vacc = u[ode_ix(c_sus, 1:k, k)] .- count_sus_to_inf

    flow_sus_vacc = sus_not_vacc[1] * flow_vacc 
    count_vacc = rand(Binomial(sus_not_vacc[1], -expm1(-flow_sus_vacc * dt)))
        

    du[ode_ix(c_sus, 1:k, k)] .+= -count_sus_to_inf
    du[ode_ix(c_inf, 1:k, k)] .+= count_sus_to_inf 

    du[ode_ix(c_count, 1:k, k)] .+= count_sus_to_inf 

    du[ode_ix(c_sus, 1, k)] -= count_vacc
    du[ode_ix(c_sus, k, k)] += count_vacc

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
    
        for i in 1:100
            ctmc_step!(du, u, model_params, 0.01)
            u = u + du
        end
    end

    return u_t, inf
end
