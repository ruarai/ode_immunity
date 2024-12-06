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



function ode_step_static(u, p)
    S = p.S
    du = MArray{S + 1}(zeros(S + 1))
    
    I = u[S + 1] 

    du[1] = -p.beta * (1 - p.p_acq[1]) * I * u[1] + 
        p.gamma * p.p_trans[1] * I +
        p.wane_transition_rate * u[2]

    for i in 2:(S - 1)
        du[i] = -p.beta * (1 - p.p_acq[i]) * I * u[i] + 
            p.gamma * p.p_trans[i] * I +
            p.wane_transition_rate * u[i + 1] -
            p.wane_transition_rate * u[i]
    end

    du[S] = -p.beta * (1 - p.p_acq[S]) * I * u[S] + 
        p.gamma * p.p_trans[S] * I -
        p.wane_transition_rate * u[S]

    du[S + 1] = sum(p.beta .* (1 .- p.p_acq[1:S]) .* I .* u[1:S])
end