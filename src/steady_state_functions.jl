

# Returns a zero vector if inf_vec implies a steady state solution
function steady_state(
    inf_vec, I_sum,
    B, M,
    omega_inv,
    k, lambda, beta
)
    ((k * lambda) / (beta * I_sum)) * (B * (inf_vec .* omega_inv)) .- inf_vec .+ M * inf_vec
end


# Gets the total population across S and I that is implied by inf_vec
function population_sum(
    inf_vec, I_sum,
    omega_inv,
    beta, gamma
)
    sus_vec = (gamma / (beta * I_sum)) * (inf_vec .* omega_inv)
    
    S_sum = sum(sus_vec)

    return S_sum + I_sum
end

# Returns a zero vector iff inf_vec is positive
function soft_inf_positive(inf_vec)
    return sqrt.(inf_vec .^ 2) - inf_vec
end

function steady_state_and_valid(
    inf_vec, 
    B, M,
    omega_inv,
    k, lambda, beta, gamma
)
    I_sum = sum(inf_vec)

    a = steady_state(inf_vec, I_sum, B, M, omega_inv, k, lambda, beta)
    b = 1 - population_sum(inf_vec, I_sum, omega_inv, beta, gamma)
    c = soft_inf_positive(inf_vec)

    a .^ 2 .+ b ^ 2 .+ c .^ 2
end



function get_steady_state(model_params)
    omega_inv = 1 ./ (1 .- model_params.p_acq)

    N = model_params.N

    u0 = SVector{N}(convert.(Double64, [1.0 / N for x in 1:N]))

    fn_solve(inf_vec, p) = steady_state_and_valid(
        inf_vec, 
        model_params.B, model_params.M, 
        omega_inv, 
        model_params.k, model_params.lambda, model_params.beta, model_params.gamma
    )
    
    probN = NonlinearProblem(fn_solve, u0)
    inf_vec = solve(
        probN, NewtonRaphson();
        abstol = 1e-30, maxiters = 4000,
        show_trace = Val(true), trace_level = NonlinearSolve.TraceAll(10)
    )

    sus_vec = (model_params.gamma / (model_params.beta * sum(inf_vec))) * (inf_vec .* omega_inv)


    return (sus_vec, inf_vec)
end