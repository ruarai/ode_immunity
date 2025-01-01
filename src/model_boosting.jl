function ode_step_boosting!(du, u, model_params, t)
    S = model_params.S

    for i in 1:(S * 3)
        du[i] = 0
    end
    flow_decay = @views u[ode_ix_boosting(c_sus, S, 2:S)] .* model_params.wane_transition_rate

    du[ode_ix_boosting(c_sus, S, 2:S)] = -flow_decay
    du[ode_ix_boosting(c_sus, S, 1:(S - 1))] .+= flow_decay

    beta_t = model_params.beta * (1 + model_params.eta * sin(2 * pi * t / 365.0 + 0.5 * pi))

    flow_sus_to_inf = @views u[ode_ix_boosting(c_sus, S, 1:S)] .*
        sum(u[ode_ix_boosting(c_inf, S, 1:S)]) .*
        beta_t .* 
        (1 .- model_params.p_acq)

    flow_inf_to_sus = @views u[ode_ix_boosting(c_inf, S, 1:S)] .* model_params.gamma


    du[ode_ix_boosting(c_sus, S, 1:S)] .+= model_params.M * flow_inf_to_sus - flow_sus_to_inf
    du[ode_ix_boosting(c_inf, S, 1:S)] .+= flow_sus_to_inf - flow_inf_to_sus

    du[ode_ix_boosting(c_count, S, 1:S)] .+= flow_sus_to_inf 
end

function ode_step_no_count_boosting!(du, u, model_params, t)
    S = model_params.S

    for i in 1:(S * 2)
        du[i] = 0
    end
    flow_decay = u[ode_ix_boosting(1, S, 2:S)] .* model_params.wane_transition_rate

    du[ode_ix_boosting(1, S, 2:S)] = -flow_decay
    du[ode_ix_boosting(1, S, 1:(S - 1))] .+= flow_decay

    beta_t = model_params.beta * (1 + model_params.eta * sin(2 * pi * t / 365.0 + 0.5 * pi))

    flow_sus_to_inf = @views u[ode_ix_boosting(c_sus, S, 1:S)] .*
        sum(u[ode_ix_boosting(c_inf, S, 1:S)]) .*
        beta_t .* 
        (1 .- model_params.p_acq)

    flow_inf_to_sus = @views u[ode_ix_boosting(c_inf, S, 1:S)] .* model_params.gamma


    du[ode_ix_boosting(c_sus, S, 1:S)] .+= model_params.M * flow_inf_to_sus - flow_sus_to_inf
    du[ode_ix_boosting(c_inf, S, 1:S)] .+= flow_sus_to_inf - flow_inf_to_sus 
end