println("Running: period over grid")

arg_ix = parse(Int, ARGS[1])
n_array = parse(Int, ARGS[2])
println("Job at array index $arg_ix of $n_array, with n_cpu = $(Threads.nthreads())")

include("../dependencies.jl")


k = 32
C = 8.0
R = 1.5
gamma = 0.25
beta = R * gamma
rho = 0.003
eta_0 = 0.0
b = 2^3
h = 8
c_jump_dist = Normal(2^6, 2^5)


model_params_0 = make_model_parameters(
    k = k, beta = beta, gamma = gamma, C = C, rho = rho,
    b = b, h = h, c_jump_dist = c_jump_dist; boosting = "none", eta = eta_0
)

ode_sparsity = ode_get_sparsity(model_params_0)


n_inf_0 = 0.0001

ϵ = 1e-6
Δt = 0.25
t = n_days_burn_in:Δt:n_days
t_daily = (n_days_burn_in + 1):n_days

# Seasonality
x_eta = 0.00:0.001:0.5
length(x_eta)

# Waning
rho_step = 0.00001
x_rho = rho_step:rho_step:0.005
length(x_rho)

x_vals = vec([(x1, x2) for x1 in x_eta, x2 in x_rho])

ix_jobs = get_jobs(arg_ix, n_array, length(x_vals))
x_vals_job = x_vals[ix_jobs]



y_period = zeros(length(x_vals_job), 3)
y_inf_summary = zeros(length(x_vals_job), 6)
y_attack_rate = zeros(length(x_vals_job))

time_start = Base.time()

Threads.@threads for i in eachindex(x_vals_job)
    println("Running job $(ix_jobs[i])")

    model_params = make_model_parameters(
        k = k, beta = beta, gamma = gamma, C = C, rho = x_vals_job[i][2],
        b = b, h = h, c_jump_dist = c_jump_dist; boosting = "none", eta = x_vals_job[i][1]
    )

    ode_solution = ode_solve(model_params, n_days, n_inf_0, ode_sparsity, saveat = Δt)

    y_inf = vec(sum(ode_solution(t)[ode_ix(c_inf, 1:model_params.S, model_params.S), :], dims = 1))
    y_inf_summary[i, 1] = minimum(y_inf)
    y_inf_summary[i, 2] = maximum(y_inf)
    y_inf_summary[i, 3] = mean(y_inf)

    y_count = vec(sum(ode_solution(t)[ode_ix(c_count, 1:model_params.S, model_params.S), :], dims = 1))
    y_inc = diff(y_count)
    y_inf_summary[i, 4] = minimum(y_inc)
    y_inf_summary[i, 5] = maximum(y_inc)
    y_inf_summary[i, 6] = mean(y_inc)


    period_mean, period_sd, period_n = get_period(ode_solution, model_params, n_days_burn_in, n_days, Δt, ϵ)

    y_period[i, :] = [period_mean period_sd period_n]
    y_attack_rate[i] = get_periodic_attack_rate(ode_solution, model_params, n_days_burn_in, n_days, period_mean)
end

time_elapsed = Base.time() - time_start

println("Completed $(length(ix_jobs)) jobs in $(round(time_elapsed, digits = 2)), ($(round(time_elapsed/length(ix_jobs), digits = 2)) seconds/job)")

x_vals_job = stack(x_vals_job)
jldsave("data_dist/period_grid/$(arg_ix).jld2"; x_vals_job, y_period, y_inf_summary, y_attack_rate)

println("Outputs saved.")