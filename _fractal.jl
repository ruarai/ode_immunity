

include("dependencies.jl")



lambdas = 0.0:0.1:0.8



for lambda in lambdas

    println(lambda)


    model_params = make_model_parameters(
        c_max = 8, k = 2,
    
        beta = 1.0, sigma = 1.0, gamma = 0.5, lambda = lambda,
    
        b = 4, m = 4,
    
        c_jump_dist = Normal(6, 0.5)
    )
    
    
    
    
    
    solver = get_steady_state(model_params)
    
    inf_sol = solver.u
    exp_sol = model_params.gamma / model_params.sigma * inf_sol
    
    I_sol = sum(inf_sol)
    
    omega_inv = 1 ./ (1 .- model_params.p_acq)
    
    sus_sol = (model_params.gamma * inf_sol .* omega_inv) / (model_params.beta * I_sol)
    
    
    u0 = zeros(Float64x4, 3 * model_params.N)
    u0_vec = vcat(sus_sol, exp_sol, inf_sol)
    
    for i in eachindex(u0)
        u0[i] = u0_vec[i]
    end
    
    using ForwardDiff
    
    du0 = copy(u0)
    J = ForwardDiff.jacobian((du, u) -> ode_step_no_count!(du, u, model_params, 0.0), du0, u0)
    J_eigvals = eigvals(J)
    
    
    
    real(J_eigvals)[end]
    scatter(((J_eigvals)))

    println(real(J_eigvals)[end])
end
