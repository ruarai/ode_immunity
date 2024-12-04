
include("dependencies.jl")


using ChaosTools
function get_max_lyapunov(model_params, n_days, n_days_burn_in)
    u0 = zeros(Float64, n_compartments * model_params.S + 2)
    u0[ode_ix(c_sus, 1, model_params.S)] = 1.0 - 0.01
    u0[ode_ix(c_inf, 1, model_params.S)] = 0.01

    u0[end - 1] = model_params.eta
    u0[end] = 0
    
    
    dyn_system = CoupledODEs(ode_step_no_count!, u0, model_params)
    return lyapunov(dyn_system, n_days - n_days_burn_in, Ttr = n_days_burn_in, Δt = 100.0, d0 = 1e-9)
end

x_eta = 0:0.01:0.5
y_lyapunov = zeros(length(x_eta))

@showprogress Threads.@threads for i_eta in eachindex(x_eta)
    model_params = make_model_parameters(
        k = baseline_k, beta = baseline_beta, gamma = baseline_gamma,
        C = baseline_C, r = 0.05,
        b = baseline_b, h = baseline_h, c_jump_dist = baseline_c_jump_dist;
        eta = x_eta[i_eta] 
    )

    y_lyapunov[i_eta] = get_max_lyapunov(model_params, n_days, n_days_burn_in)
end


plot(x_eta, y_lyapunov)
x_eta[y_lyapunov .> 0.0002]



model_params = make_model_parameters(
    k = baseline_k, beta = baseline_beta, gamma = baseline_gamma,
    C = baseline_C, r = 0.05,
    b = baseline_b, h = baseline_h, c_jump_dist = baseline_c_jump_dist;
    eta = 0.5
)



u0 = zeros(Float64, n_compartments * model_params.S + 2)
u0[ode_ix(c_sus, 1, model_params.S)] = 1.0 - 0.01
u0[ode_ix(c_inf, 1, model_params.S)] = 0.01

u0[end - 1] = model_params.eta
u0[end] = 0


dyn_system = CoupledODEs(ode_step_no_count!, u0, model_params)

s = gali(dyn_system, 40000, 3; u0 = u0)[2][end]

predictability(dyn_system)


return lyapunov(dyn_system, n_days - n_days_burn_in, Ttr = n_days_burn_in, Δt = 100.0, d0 = 1e-9)