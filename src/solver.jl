

function ode_ix_boosting(c, S, i)
    return (c .- 1) .* S .+ i
end

function ode_ix_sus(i)
    return i
end

function ode_ix_inf(S)
    return S + 1
end

function ode_ix_count(S)
    return S + 2
end


function ode_get_sparsity(
    model_params
)
    u0 = zeros(Float64, model_params.S + 2)
    du0 = copy(u0)
    return Symbolics.jacobian_sparsity((du, u) -> ode_step!(du, u, model_params, 0.0), du0, u0)
end


function ode_solve(
    model_params,
    n_days,
    n_inf_0,
    ode_sparsity;
    saveat_step = 1,
    datatype = Float64,
    n_days_burn_in = 0.0
)
    ode_step_fn = ODEFunction(ode_step!; jac_prototype = float.(ode_sparsity))

    u0 = zeros(datatype, model_params.S + 2)
    u0[ode_ix_sus(1)] = 1.0 - n_inf_0
    u0[ode_ix_inf(model_params.S)] = n_inf_0

    tspan = (0.0, n_days)
    prob = ODEProblem{true, SciMLBase.FullSpecialize}(ode_step_fn, u0, tspan, model_params)

    vec_saveat = n_days_burn_in:saveat_step:n_days

    return DifferentialEquations.solve(
        prob, Rodas4P(), 
        dtmax = 8.0,
        abstol = 1e-16,
        reltol = 1e-5,
        saveat = vec_saveat
    );
end