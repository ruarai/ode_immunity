
# Return a zero vector if u is positive
function soft_inf_positive(u)
    return sqrt.(u .^ 2) - u
end

function steady_state_and_valid(u, model_params)
    du = zeros(eltype(u), length(u))
    ode_step_no_count!(du, u, model_params, 0)

    a = du # Is total change zero?
    b = 1 - sum(u) # Is population total size 1?
    c = soft_inf_positive(u) # Is each compartment positive valued?

    a .^ 2 .+ b ^ 2 .+ c .^ 2
end


function steady_state_and_valid_boosting(u, model_params_boosting)
    du = zeros(eltype(u), length(u))
    ode_step_no_count_boosting!(du, u, model_params_boosting, 0)

    a = du # Is total change zero?
    b = 1 - sum(u) # Is population total size 1?
    c = soft_inf_positive(u) # Is each compartment positive valued?

    a .^ 2 .+ b ^ 2 .+ c .^ 2
end


function get_steady_state(model_params, verbose = false)
    n_compartments = model_params.S + 1

    # 2.0 works for whatever reason.
    u0 = convert.(Float64, [2.0 / n_compartments for x in 1:n_compartments])
    
    fn_solve(u, p) = steady_state_and_valid(u, model_params)
    nonlinear_prob = NonlinearProblem(fn_solve, u0)
    
    try
        u_sol = solve(
            nonlinear_prob, NonlinearSolve.NewtonRaphson();
            abstol = 1e-30, maxiters = 8000,
            show_trace = Val(verbose), trace_level = NonlinearSolve.TraceAll(10)
        )
        
        return u_sol
    catch e
        println("Failed to find steady state with exception $e")

        return zeros(n_compartments)
    end
end


function get_steady_state_boosting(model_params_boosting, verbose = false)
    n_compartments = model_params_boosting.S * 2

    # 2.0 works for whatever reason.
    u0 = convert.(Float64, [2.0 / n_compartments for x in 1:n_compartments])
    
    fn_solve(u, p) = steady_state_and_valid_boosting(u, model_params_boosting)
    nonlinear_prob = NonlinearProblem(fn_solve, u0)
    
    try
        u_sol = solve(
            nonlinear_prob, NonlinearSolve.NewtonRaphson();
            abstol = 1e-30, maxiters = 8000,
            show_trace = Val(verbose), trace_level = NonlinearSolve.TraceAll(10)
        )
        
        return u_sol
    catch e
        println("Failed to find steady state with exception $e")

        return zeros(n_compartments)
    end
end