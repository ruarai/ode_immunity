function ode_step!(du, u, model_params, t)
    S = model_params.S

    for i in 1:(S + 2)
        du[i] = 0
    end
    flow_decay = @views u[ode_ix_sus(2:S)] .* model_params.wane_transition_rate

    du[ode_ix_sus(2:S)] = -flow_decay
    du[ode_ix_sus(1:(S - 1))] .+= flow_decay

    beta_t = model_params.beta * (1 + model_params.eta * sin(2 * pi * t / 365.0 + 0.5 * pi))

    flow_sus_to_inf = @views u[ode_ix_sus(1:S)] .*
        sum(u[ode_ix_inf(S)]) .*
        beta_t .* 
        (1 .- model_params.p_acq)

    inf_incidence = sum(flow_sus_to_inf)

    flow_inf_to_sus = u[ode_ix_inf(S)] * model_params.gamma



    du[ode_ix_sus(1:S)] .+= model_params.p_trans * flow_inf_to_sus - flow_sus_to_inf
    du[ode_ix_inf(S)] += inf_incidence - flow_inf_to_sus

    du[ode_ix_count(S)] += inf_incidence 
end

function ode_step_no_count!(du, u, model_params, t)
    S = model_params.S

    for i in 1:(S + 1)
        du[i] = 0
    end
    flow_decay = @views u[ode_ix_sus(2:S)] .* model_params.wane_transition_rate

    du[ode_ix_sus(2:S)] = -flow_decay
    du[ode_ix_sus(1:(S - 1))] .+= flow_decay

    beta_t = model_params.beta * (1 + model_params.eta * sin(2 * pi * t / 365.0 + 0.5 * pi))

    flow_sus_to_inf = @views u[ode_ix_sus(1:S)] .*
        sum(u[ode_ix_inf(S)]) .*
        beta_t .* 
        (1 .- model_params.p_acq)

    inf_incidence = sum(flow_sus_to_inf)

    flow_inf_to_sus = u[ode_ix_inf(S)] * model_params.gamma


    du[ode_ix_sus(1:S)] .+= model_params.p_trans * flow_inf_to_sus - flow_sus_to_inf
    du[ode_ix_inf(S)] += inf_incidence - flow_inf_to_sus
end