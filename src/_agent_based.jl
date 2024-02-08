


include("dependencies.jl")

using JLD2

k = 16
R = 1.5
gamma = 0.25
beta = R * gamma
b = 0.25
lambda = 0.008
m = 40
c_jump_dist = Normal(0.8, 0.05)

# Temporary model_parameters so sparsity can be generated
model_params = make_model_parameters(
    k = k, beta = beta, gamma = gamma, lambda = lambda,
    b = b, m = m, c_jump_dist = c_jump_dist; boosting = true
)



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
                p_decay = -expm1(- model_params.lambda * k * dt)

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




n_days = 1000

n_pop_size = 100000
n_track = 100


track_state, I_t = simulate_agent_based(n_days, n_track, n_pop_size, model_params, 2)

track_strata = (track_state .- 1) .% k .+ 1
track_compartment = (track_state .> ode_ix(c_sus, k, k)) .* 100
plot(track_strata[:,10:20], lc = track_compartment[:,10:20], legend = false)


plot(1:n_days, I_t)


jldsave("data/paper/agent_based.jld2"; track_state, I_t)