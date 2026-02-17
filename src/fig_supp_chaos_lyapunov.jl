
include("dependencies.jl")

model_params_base = make_model_parameters(
    k = baseline_k, beta = baseline_beta, gamma = baseline_gamma,
    a = baseline_a, r = 0.018,
    b = baseline_b, h = baseline_h, c_jump_dist = baseline_c_jump_dist;
)

n_days_burnin = 365 * 1000
n_days = n_days_burnin + 365 * 250
seq_t = collect(n_days_burnin:n_days)

ode_step_fn = ODEFunction(ode_step!)

u0 = zeros(Float64, model_params_base.S + 2)
u0[ode_ix_sus(1)] = 1.0 - n_inf_0
u0[ode_ix_inf(model_params_base.S)] = n_inf_0



x_eta = collect(0.00:0.001:0.5)
x_r = collect(0.01:0.00006:0.025)

x_vals = vec([(eta = x1, r = x2) for x1 in x_eta, x2 in x_r])


y_lyapunov = zeros(length(x_vals))

@showprogress Threads.@threads for i in eachindex(x_vals)
    model_params = make_model_parameters(
        k = baseline_k, beta = baseline_beta, gamma = baseline_gamma,
        a = baseline_a, r = x_vals[i].r,
        b = baseline_b, h = baseline_h, c_jump_dist = baseline_c_jump_dist;

        eta = x_vals[i].eta
    )

    prob = ContinuousDynamicalSystem(ode_step!, u0, model_params; t0 = 0.0)
    y_lyapunov[i] = lyapunov(prob, n_days, Ttr = n_days_burnin, Δt = 100)
end


y = reshape(y_lyapunov, length(x_eta), length(x_r))

heatmap(x_eta, x_r, y')

x_eta_full = [x_vals[i].eta for i in eachindex(x_vals)]
x_r_full = [x_vals[i].r for i in eachindex(x_vals)]

jldsave("data/chaos_lyapunov.jld2"; x_eta_full, x_r_full, y_lyapunov)
