
include("dependencies.jl")


using BifurcationKit, Accessors


function ode_step_bifur!(du, u, model_params)
    ode_step_minimal!(du, u, model_params, 0.0)
    return du
end

function rec_sol(u, p; params...)
    return 1 - sum(u)
end

prob = BifurcationProblem(
    ode_step_bifur!, x0.u[1:33], model_params, (@optic _.wane_transition_rate);

    record_from_solution = rec_sol
)

opts = ContinuationPar(
    p_min = 0.0, p_max = 0.18,
    ds = 1e-5, dsmax = 1e-5, dsmin = 1e-7,           
    detect_bifurcation = 3,    
    max_steps = 20000,
    nev = 33
)

br = continuation(prob, PALC(), opts,)


# Lower bifurcation
get_normal_form(br, 1)

# Upper bifurction
get_normal_form(br, 2)

plot(br, vars = (:param, :x))