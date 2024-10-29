

function ode_ix(c, n, k)
    return (c .- 1) .* k .+ n
end


function ode_step!(du, u, model_params, t)
    k = model_params.S

    for i in 1:(k * n_compartments)
        du[i] = 0
    end

    flow_decay = @views u[ode_ix(1, 2:k, k)] .* model_params.wane_transition_rate

    du[ode_ix(1, 2:k, k)] = -flow_decay
    du[ode_ix(1, 1:(k - 1), k)] .+= flow_decay

    beta_t = model_params.beta * (1 + model_params.eta * sin(2 * pi * t / 365.0))

    flow_sus_to_inf = @views u[ode_ix(c_sus, 1:k, k)] .*
        sum(u[ode_ix(c_inf, 1:k, k)]) .*
        beta_t .* 
        (1 .- model_params.p_acq)

    flow_inf_to_sus = @views u[ode_ix(c_inf, 1:k, k)] .* model_params.gamma


    du[ode_ix(c_sus, 1:k, k)] .+= model_params.M * flow_inf_to_sus - flow_sus_to_inf
    du[ode_ix(c_inf, 1:k, k)] .+= flow_sus_to_inf - flow_inf_to_sus

    du[ode_ix(c_count, 1:k, k)] .+= flow_sus_to_inf 
end


function ode_step_no_count!(du, u, model_params, t)
    k = model_params.S

    for i in 1:(k * 3)
        du[i] = 0
    end

    flow_decay = u[ode_ix(1, 2:k, k)] .* model_params.wane_transition_rate

    du[ode_ix(1, 2:k, k)] = -flow_decay
    du[ode_ix(1, 1:(k - 1), k)] .+= flow_decay

    flow_sus_to_inf = u[ode_ix(c_sus, 1:k, k)] .*
        sum(u[ode_ix(c_inf, 1:k, k)]) .*
        model_params.beta .* 
        (1 .- model_params.p_acq)
    
    flow_inf_to_sus = u[ode_ix(c_inf, 1:k, k)] .* model_params.gamma


    du[ode_ix(c_sus, 1:k, k)] .+= model_params.M * flow_inf_to_sus - flow_sus_to_inf
    du[ode_ix(c_inf, 1:k, k)] .+= flow_sus_to_inf - flow_inf_to_sus
    
end


function ode_get_sparsity(
    model_params
)
    u0 = zeros(Float64, n_compartments * model_params.S)
    du0 = copy(u0)
    return Symbolics.jacobian_sparsity((du, u) -> ode_step!(du, u, model_params, 0.0), du0, u0)
end


function ode_solve(
    model_params,
    n_days,
    n_inf_0,
    ode_sparsity;
    saveat = 1,
    dt = 0.01,
    datatype = Float64
)
    ode_step_fn = ODEFunction(ode_step!; jac_prototype = float.(ode_sparsity))

    u0 = zeros(datatype, n_compartments * model_params.S)
    u0[ode_ix(c_sus, 1, model_params.S)] = 1.0 - n_inf_0
    u0[ode_ix(c_inf, 1, model_params.S)] = n_inf_0

    tspan = (0.0, 1.0 * n_days)
    prob = ODEProblem{true, SciMLBase.FullSpecialize}(ode_step_fn, u0, tspan, model_params)

    return DifferentialEquations.solve(prob, Euler(), dt = dt, saveat = saveat);
end