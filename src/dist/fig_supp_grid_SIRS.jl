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

sigma_step = 0.00002
x_sigma = sigma_step:sigma_step:0.01
length(x_sigma)

length(x_eta) * length(x_sigma)

x_vals = vec([(eta = x1, sigma = x2) for x1 in x_eta, x2 in x_sigma])

ix_jobs = get_jobs(arg_ix, n_array, length(x_vals))
x_vals_job = x_vals[ix_jobs]

y_period = zeros(length(x_vals_job), 3)
y_inf_summary = zeros(length(x_vals_job), 11)

function sir_ode!(du, u, p, t)
    S, I, R = u

    beta_t = p.beta * (1 + p.eta * cos(2π * t / 365.0))

    du[1] = -beta_t * S * I + p.sigma * R
    du[2] = beta_t * S * I - p.gamma * I
    du[3] = p.gamma * I - p.sigma * R
end
u0 = [0.99, 0.01, 0.0]

time_start = Base.time()


@showprogress Threads.@threads for i in eachindex(x_vals_job)
    println("Running job $(ix_jobs[i])")
    p = (beta = baseline_beta, gamma = baseline_gamma, sigma = x_vals_job[i].sigma, eta = x_vals_job[i].eta)
    tspan = (0.0, n_days) 

    prob = ODEProblem(sir_ode!, u0, tspan, p)

    sol = solve(
        prob,
        Rodas5P(),

        dtmax = 1.0,
        reltol = 1e-10, abstol = 1e-10,

        saveat = n_days_burn_in:0.25:n_days
    )

    y_period[i, :] .= get_period(sol, n_days_burn_in, n_days, 0.25, periodic_ϵ, 1:3)

    inf = [sol(t)[2] for t in n_days_burn_in:n_days]

    y_inf_summary[i, 1] = minimum(inf)
    y_inf_summary[i, 2] = maximum(inf)
    y_inf_summary[i, 3] = mean(inf)
    y_inf_summary[i, 4] = testchaos01(NaNMath.log10.(inf[1:80:end]))

    # prob_lya = ContinuousDynamicalSystem(sir_ode!, u0, p; t0 = 0.0, diffeq = (alg = Rodas5P(), abstol = 1e-10, reltol = 1e-10))
    # y_inf_summary[i, 11] = lyapunov(prob_lya, n_days, Ttr = n_days_burn_in, Δt = 100)
end

time_elapsed = Base.time() - time_start

println("Completed $(length(ix_jobs)) jobs in $(round(time_elapsed, digits = 2)), ($(round(time_elapsed/length(ix_jobs), digits = 2)) seconds/job)")

x_vals_job = stack(x_vals_job)
mkpath("data_dist/period_grid_SIRS")
jldsave("data_dist/period_grid_SIRS/$(arg_ix).jld2"; x_vals_job, y_period, y_inf_summary)

println("Outputs saved.")