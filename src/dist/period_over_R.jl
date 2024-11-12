println("Running: period over grid")

include("../dependencies.jl")

arg_ix = parse(Int, ARGS[1])
n_array = parse(Int, ARGS[2])
println("Job at array index $arg_ix of $n_array, with n_cpu = $(Threads.nthreads())")


k = 32
C = 8.0
R_0 = 1.5
gamma = 0.25
beta_0 = R_0 * gamma
rho_0 = 0.003
b = 2^3
h = 8
c_jump_dist = Normal(2^6, 2^5)


model_params_0 = make_model_parameters(
    k = k, beta = beta, gamma = gamma, C = C, rho = rho_0,
    b = b, h = h, c_jump_dist = c_jump_dist; boosting = "none"
)

ode_sparsity = ode_get_sparsity(model_params_0)


n_inf_0 = 0.0001
n_days_burn_in = 30000
n_days = (n_days_burn_in + 365 * 20 * 2)

ϵ = 1e-6
Δt = 0.25
t = n_days_burn_in:Δt:n_days
t_daily = (n_days_burn_in + 1):n_days

# R0
x_R0 = exp.(0:0.01:0.7)
length(x_R0)

# Waning
rho_step = 0.0005
x_rho = rho_step:rho_step:0.005
length(x_rho)

x_vals = vec([(R0 = x1, rho = x2) for x1 in x_R0, x2 in x_rho])

ix_jobs = get_jobs(arg_ix, n_array, length(x_vals))
x_vals_job = x_vals[ix_jobs]



y_period = zeros(length(x_vals_job), 3)
y_inf_maxima = zeros(length(x_vals_job), 3)
y_attack_rate = zeros(length(x_vals_job))

time_start = Base.time()

Threads.@threads for i in eachindex(x_vals_job)
    println("Running job $(ix_jobs[i])")

    beta_i = x_vals_job[i].R0 * gamma
    model_params = make_model_parameters(
        k = k, beta = beta_i, gamma = gamma, C = C, rho = x_vals_job[i].rho,
        b = b, h = h, c_jump_dist = c_jump_dist; boosting = "none"
    )

    ode_solution = ode_solve(model_params, n_days, n_inf_0, ode_sparsity, saveat = Δt)

    y = ode_solution(t)[1:(model_params.S * 2), :]'

    y_inf = ode_solution(t)[(model_params.S + 1):(model_params.S * 2), :]'
    y_inf_sum = sum(y_inf, dims = 2)
    y_inf_maxima[i, 1] = minimum(y_inf_sum)
    y_inf_maxima[i, 2] = maximum(y_inf_sum)
    y_inf_maxima[i, 3] = mean(y_inf_sum)

    period_mean, period_sd, period_n = get_period(ode_solution, model_params, n_days_burn_in, n_days, Δt, ϵ)

    y_period[i, :] = [period_mean period_sd period_n]
    y_attack_rate[i] = get_periodic_attack_rate(ode_solution, model_params, n_days_burn_in, n_days, period_mean)
end


plot(y_period)

heatmap(x_rho, x_R0, reshape(y_period[:,1], (length(x_R0), length(x_rho))))