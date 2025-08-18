

n_days_burn_in = 1000 * 365
n_days = n_days_burn_in + 250 * 365
t_post_burn_in = collect(n_days_burn_in:n_days)

x_eta = collect(0.2:0.001:0.5)
length(x_eta)

y_inf = zeros(length(x_eta), length(t_post_burn_in))
y_inc = zeros(length(x_eta), length(t_post_burn_in))

@showprogress Threads.@threads for i in eachindex(x_eta)
    model_params = make_model_parameters(
        k = baseline_k, beta = baseline_beta, gamma = baseline_gamma,
        a = baseline_a, r = 0.018,
        b = baseline_b, h = baseline_h, c_jump_dist = baseline_c_jump_dist;
        eta = x_eta[i]
    )

    ode_solution = ode_solve(model_params, n_days, n_inf_0, n_days_burn_in = n_days_burn_in)

    sus, inf, inc = get_results(ode_solution, t_post_burn_in, model_params)

    y_inf[i, :] = inf
    y_inc[i, :] = inc
end



scatter(fill(1, 249), [y_inc[1, 365 * t + 180] for t in 1:249], legend = false, ms = 0.1)

for i in eachindex(x_eta)
    scatter!(fill(i, 249), [y_inc[i, 365 * t + 180] for t in 1:249], legend = false, ms = 0.1)
end

plot!(ylims = (0, 0.001))

downsample_rate = vcat(1, collect(10:10:200))

chaos_result = zeros(length(x_eta), length(downsample_rate))

@showprogress Threads.@threads for i in eachindex(x_eta)
    for ix_d in eachindex(downsample_rate)
        chaos_result[i, ix_d] = testchaos01(NaNMath.log10.(y_inc[i, 2:downsample_rate[ix_d]:end]))
    end
end

heatmap(x_eta, downsample_rate, chaos_result')


jldsave("data/chaos_downsample.jld2"; x_eta, t_post_burn_in, downsample_rate, y_inc, chaos_result)
