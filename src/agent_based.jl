

function step_abm!(population_state, n_pop_size, model_params, dt)
    k = model_params.k

    I = sum(population_state .> ode_ix(c_sus, k, k))

    for i in 1:n_pop_size
        state = population_state[i]
        new_state = state

        infected = state > ode_ix(c_sus, k, k)
        strata = (state - 1) % k + 1

        r = rand()

        if infected
            p_recover = -expm1(- model_params.gamma * dt)
            
            if r < p_recover
                new_strata = wsample(1:k, model_params.M[:,strata])
                new_state = ode_ix(c_sus, new_strata, k)
            end
        else
            # Susceptible

            p_infected = -expm1(- model_params.beta * (1 - model_params.p_acq[strata]) * I * dt / n_pop_size)

            if strata > 1
                p_decay = -expm1(- model_params.wane_transition_rate * dt)

                if r < p_decay
                    new_state = ode_ix(c_sus, strata - 1, k)
                elseif r < p_decay + p_infected 
                    new_state = ode_ix(c_inf, strata, k)
                end
            else
                if r < p_infected 
                    new_state = ode_ix(c_inf, strata, k)
                end
            end
        end

        population_state[i] = new_state
    end

end

function simulate_agent_based(n_days, n_track, n_pop_size, model_params, seed)

    population_state = zeros(Int8, n_pop_size)

    population_state[1:n_pop_size] .= ode_ix(c_sus, 1, k)
    population_state[1:10] .= ode_ix(c_inf, 1, k)
    
    track_state = zeros(n_days, n_track)
    I_t = zeros(n_days)

    Random.seed!(seed)
    
    for d in 1:n_days
        println(d)
        for i in 1:n_track
            track_state[d, i] = population_state[i]
        end
        I_t[d] = sum(population_state .> ode_ix(c_sus, model_params.k, model_params.k))
    
        for t in 1:100
            step_abm!(population_state, n_pop_size, model_params, 0.01)
        end
     
    end

    return track_state, I_t
end
