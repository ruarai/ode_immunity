

function ode_ix(c, n, N)
    return (c .- 1) .* N .+ n
end


function ode_step!(du, u, model_params, t)
    N = model_params.N

    for i in 1:(N * n_compartments)
        du[i] = 0
    end

    flow_decay = @views u[ode_ix(1, 2:N, N)] .* model_params.lambda .* model_params.k

    du[ode_ix(1, 2:N, N)] = -flow_decay
    du[ode_ix(1, 1:(N - 1), N)] .+= flow_decay

    flow_sus_to_inf = @views u[ode_ix(c_sus, 1:N, N)] .*
        sum(u[ode_ix(c_inf, 1:N, N)]) .*
        model_params.beta .* 
        (1 .- model_params.p_acq)

    flow_inf_to_sus = @views u[ode_ix(c_inf, 1:N, N)] .* model_params.gamma


    du[ode_ix(c_sus, 1:N, N)] .+= model_params.M * flow_inf_to_sus - flow_sus_to_inf
    du[ode_ix(c_inf, 1:N, N)] .+= flow_sus_to_inf - flow_inf_to_sus

    du[ode_ix(c_count, 1:N, N)] .+= flow_sus_to_inf 
end


function ode_step_no_count!(du, u, model_params, t)
    N = model_params.N

    for i in 1:(N * 3)
        du[i] = 0
    end

    flow_decay = u[ode_ix(1, 2:N, N)] .* model_params.lambda .* model_params.k

    du[ode_ix(1, 2:N, N)] = -flow_decay
    du[ode_ix(1, 1:(N - 1), N)] .+= flow_decay

    flow_sus_to_inf = u[ode_ix(c_sus, 1:N, N)] .*
        sum(u[ode_ix(c_inf, 1:N, N)]) .*
        model_params.beta .* 
        (1 .- model_params.p_acq)
    
    flow_inf_to_sus = u[ode_ix(c_inf, 1:N, N)] .* model_params.gamma


    du[ode_ix(c_sus, 1:N, N)] .+= model_params.M * flow_inf_to_sus - flow_sus_to_inf
    du[ode_ix(c_inf, 1:N, N)] .+= flow_sus_to_inf - flow_inf_to_sus
    
end


function ode_get_sparsity(
    model_params
)
    u0 = zeros(Float64, n_compartments * model_params.N)
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

    u0 = zeros(datatype, n_compartments * model_params.N)
    u0[ode_ix(c_sus, 1, model_params.N)] = 1.0 - n_inf_0
    u0[ode_ix(c_inf, 1, model_params.N)] = n_inf_0

    tspan = (0.0, 1.0 * n_days)
    prob = ODEProblem{true, SciMLBase.FullSpecialize}(ode_step_fn, u0, tspan, model_params)

    return DifferentialEquations.solve(prob, Euler(), dt = dt, saveat = saveat);
end