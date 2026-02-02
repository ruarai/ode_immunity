function ode_step!(du, u, model_params, t)
    S = model_params.S
    fill!(du, 0) 

    I = u[S + 1]

    # Waning transitions between susceptibility compartments
    @inbounds @simd for ix_sus in 2:S
        decay = u[ix_sus] * model_params.wane_transition_rate
        du[ix_sus] -= decay
        du[ix_sus - 1] += decay
    end

    # Seasonal forcing term
    beta_t_I = model_params.beta * (1 + model_params.eta * cos(2π * t / 365.0)) * (I + model_params.importation_rate)
    gamma_I = model_params.gamma * I

    # Infection flow from susceptible to infected
    p_acq = model_params.p_acq
    inf_incidence = 0.0
    @inbounds @simd for ix_sus in 1:S
        strata_incidence = u[ix_sus] * beta_t_I * (1 - p_acq[ix_sus])
        inf_incidence += strata_incidence

        du[ix_sus] -= strata_incidence
        du[ix_sus] += model_params.p_trans[ix_sus] * gamma_I
    end


    # Update infected and count compartments
    du[S + 1] += inf_incidence - gamma_I
    du[S + 2] += inf_incidence

    return nothing
end

function ode_step_no_count!(du, u, model_params, t)
    S = model_params.S
    fill!(du, 0)

    I = u[S + 1]

    # Waning transitions between susceptibility compartments
    @inbounds @simd for ix_sus in 2:S
        decay = u[ix_sus] * model_params.wane_transition_rate
        du[ix_sus] -= decay
        du[ix_sus - 1] += decay
    end

    # Seasonal forcing term
    beta_t_I = model_params.beta * (1 + model_params.eta * cos(2π * t / 365.0)) * (I + model_params.importation_rate)
    gamma_I = model_params.gamma * I

    # Infection flow from susceptible to infected
    p_acq = model_params.p_acq
    inf_incidence = 0.0
    @inbounds @simd for ix_sus in 1:S
        strata_incidence = u[ix_sus] * beta_t_I * (1 - p_acq[ix_sus])
        inf_incidence += strata_incidence

        du[ix_sus] -= strata_incidence
        du[ix_sus] += model_params.p_trans[ix_sus] * gamma_I
    end


    # Update infected compartment
    du[S + 1] += inf_incidence - gamma_I

    return nothing
end



function ode_step_minimal!(du, u, model_params, t)
    S = model_params.S
    fill!(du, 0)

    I = 1 - sum(u[1:S])

    # Waning transitions between susceptibility compartments
    @inbounds @simd for ix_sus in 2:S
        decay = u[ix_sus] * model_params.wane_transition_rate
        du[ix_sus] -= decay
        du[ix_sus - 1] += decay
    end

    # Seasonal forcing term
    beta_t_I = model_params.beta * (1 + model_params.eta * cos(2π * t / 365.0)) * (I + model_params.importation_rate)
    gamma_I = model_params.gamma * I

    # Infection flow from susceptible to infected
    p_acq = model_params.p_acq
    inf_incidence = 0.0
    @inbounds @simd for ix_sus in 1:S
        strata_incidence = u[ix_sus] * beta_t_I * (1 - p_acq[ix_sus])
        inf_incidence += strata_incidence

        du[ix_sus] -= strata_incidence
        du[ix_sus] += model_params.p_trans[ix_sus] * gamma_I
    end

    return nothing
end