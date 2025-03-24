
# Index functions without boosting
function ode_ix_sus(i)
    return i
end

function ode_ix_inf(S)
    return S + 1
end

function ode_ix_count(S)
    return S + 2
end

# Index function for with boosting (i.e. with I stratification)
function ode_ix_boosting(c, S, i)
    return (c .- 1) .* S .+ i
end

function ode_solve(
    model_params,
    n_days,
    n_inf_0;
    saveat_step = 1,
    datatype = Float64,
    n_days_burn_in = 0.0,
    u0 = nothing
)
    ode_step_fn = ODEFunction(ode_step!)

    if isnothing(u0)
        u0 = zeros(datatype, model_params.S + 2)
        u0[ode_ix_sus(1)] = 1.0 - n_inf_0
        u0[ode_ix_inf(model_params.S)] = n_inf_0
    else
        u0 = vcat(u0, 0.0)
    end

    tspan = (0.0, n_days)
    prob = ODEProblem{true, SciMLBase.FullSpecialize}(ode_step_fn, u0, tspan, model_params)

    vec_saveat = n_days_burn_in:saveat_step:n_days

    return DifferentialEquations.solve(
        prob, Rodas5P(), dt = 0.01,

        dtmax = 1.0,
        reltol = 1e-10, abstol = 1e-10,

        saveat = vec_saveat
    );
end


function ode_solve_boosting(
    model_params,
    n_days,
    n_inf_0;
    saveat_step = 1,
    datatype = Float64,
    n_days_burn_in = 0.0
)
    ode_step_fn = ODEFunction(ode_step_boosting!)

    u0 = zeros(datatype, model_params.S * 2 + 1)
    u0[ode_ix_boosting(c_sus, model_params.S, 1)] = 1.0 - n_inf_0
    u0[ode_ix_boosting(c_inf, model_params.S, 1)] = n_inf_0

    tspan = (0.0, n_days)
    prob = ODEProblem{true, SciMLBase.FullSpecialize}(ode_step_fn, u0, tspan, model_params)

    vec_saveat = n_days_burn_in:saveat_step:n_days

    return DifferentialEquations.solve(
        prob, Rodas5P(), dt = 0.01,

        dtmax = 1.0,
        reltol = 1e-10, abstol = 1e-10,
        
        saveat = vec_saveat
    );
end