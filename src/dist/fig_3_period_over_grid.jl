println("Running: period over grid")

arg_ix = parse(Int, ARGS[1])
n_array = parse(Int, ARGS[2])
println("Job at array index $arg_ix of $n_array, with n_cpu = $(Threads.nthreads())")

include("../dependencies.jl")

n_days_burn_in = 1000 * 365
n_days = n_days_burn_in + 250 * 365

periodic_Δt = 0.25
t_post_burn_in = n_days_burn_in:n_days

# Seasonality
x_eta = 0.00:0.001:0.5
length(x_eta)

# Waning
r_step = 0.00006
x_r = r_step:r_step:0.03
length(x_r)

x_vals = vec([(eta = x1, r = x2) for x1 in x_eta, x2 in x_r])

ix_jobs = get_jobs(arg_ix, n_array, length(x_vals))
x_vals_job = x_vals[ix_jobs]

y_period = zeros(length(x_vals_job), 3)
y_inf_summary = zeros(length(x_vals_job), 11)
y_seasonality = zeros(length(x_vals_job), 2)

time_start = Base.time()

Threads.@threads for i in eachindex(x_vals_job)
    println("Running job $(ix_jobs[i])")

    model_params = make_model_parameters(
        k = baseline_k, beta = baseline_beta, gamma = baseline_gamma,
        a = baseline_a, r = x_vals_job[i].r,
        b = baseline_b, h = baseline_h, c_jump_dist = baseline_c_jump_dist;
        eta = x_vals_job[i].eta
    )

    ode_solution = ode_solve(
        model_params, n_days, n_inf_0,
        saveat_step = periodic_Δt, n_days_burn_in = n_days_burn_in
    )

    inf = get_inf(ode_solution, t_post_burn_in, model_params)
    inc = get_inc(ode_solution, t_post_burn_in, model_params; maintain_length = false)

    y_inf_summary[i, 1] = minimum(inf)
    y_inf_summary[i, 2] = maximum(inf)
    y_inf_summary[i, 3] = mean(inf)
    y_inf_summary[i, 4] = testchaos01(NaNMath.log10.(inf[1:80:end]))

    y_inf_summary[i, 5] = minimum(inc)
    y_inf_summary[i, 6] = maximum(inc)
    y_inf_summary[i, 7] = mean(inc)
    y_inf_summary[i, 8] = testchaos01(NaNMath.log10.(inc[1:80:end]))
    y_inf_summary[i, 9] = Integer(ode_solution.retcode)
    y_inf_summary[i, 10] = get_peak_density(inc, t_post_burn_in)
    y_inf_summary[i, 11] = entropy(inc)

    period_mean, period_sd, period_n = get_period(
        ode_solution, model_params, n_days_burn_in, n_days, 
        periodic_Δt, periodic_ϵ
    )

    y_period[i, :] = [period_mean period_sd period_n]

    seasonality_x, seasonality_y = get_seasonality_coordinates(ode_solution, t_post_burn_in, model_params)

    y_seasonality[i, :] = [seasonality_x, seasonality_y]
end

time_elapsed = Base.time() - time_start

println("Completed $(length(ix_jobs)) jobs in $(round(time_elapsed, digits = 2)), ($(round(time_elapsed/length(ix_jobs), digits = 2)) seconds/job)")

x_vals_job = stack(x_vals_job)
jldsave("data_dist/period_grid/$(arg_ix).jld2"; x_vals_job, y_period, y_inf_summary, y_seasonality)

println("Outputs saved.")