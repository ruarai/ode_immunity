function ode_step_boosting!(du, u, model_params, t)
    S = model_params.S

    fill!(du, 0)
    
    # Waning transitions between susceptibility compartments
    @inbounds @simd for ix_sus in 2:S
        decay = u[ix_sus] * model_params.wane_transition_rate
        du[ix_sus] -= decay
        du[ix_sus - 1] += decay
    end

    # Seasonal forcing term
    beta_t = model_params.beta * (1 + model_params.eta * cos(2π * t / 365.0))
    I = sum(u[(S + 1):(2 * S)])

    # Infection flow from susceptible to infected
    inf_incidence = 0.0
    @inbounds @simd for ix_sus in 1:S
        ix_inf = ix_sus + S

        strata_incidence = u[ix_sus] * beta_t * I * (1 - model_params.p_acq[ix_sus])
        inf_incidence += strata_incidence

        du[ix_sus] -= strata_incidence
        du[ix_inf] += strata_incidence - model_params.gamma * u[ix_inf]

        for i in 1:S
            du[ix_sus] += u[i + S] * model_params.M[ix_sus, i] * model_params.gamma 
        end
    end

    du[ode_ix_boosting(c_count, S, 1:S)] .+= inf_incidence 
end

function ode_step_no_count_boosting!(du, u, model_params, t)
    S = model_params.S
    fill!(du, 0)
    
    # Waning transitions between susceptibility compartments
    @inbounds @simd for ix_sus in 2:S
        decay = u[ix_sus] * model_params.wane_transition_rate
        du[ix_sus] -= decay
        du[ix_sus - 1] += decay
    end

    # Seasonal forcing term
    beta_t = model_params.beta * (1 + model_params.eta * cos(2π * t / 365.0))
    I = sum(u[(S + 1):(2 * S)])

    # Infection flow from susceptible to infected
    inf_incidence = 0.0
    @inbounds @simd for ix_sus in 1:S
        ix_inf = ix_sus + S

        strata_incidence = u[ix_sus] * beta_t * I * (1 - model_params.p_acq[ix_sus])
        inf_incidence += strata_incidence

        du[ix_sus] -= strata_incidence
        du[ix_inf] += strata_incidence - model_params.gamma * u[ix_inf]

        for i in 1:S
            du[ix_sus] += u[i + S] * model_params.M[ix_sus, i] * model_params.gamma 
        end
    end
end