println("Running: period over grid")

include("../dependencies.jl")

arg_ix = parse(Int, ARGS[1])
n_array = parse(Int, ARGS[2])
println("Job at array index $arg_ix of $n_array, with n_cpu = $(Threads.nthreads())")


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
n_days_burn_in = 30000
n_days = (n_days_burn_in + 365 * 100 * 2)

ϵ = 1e-6
Δt = 0.25
t = n_days_burn_in:Δt:n_days
t_daily = (n_days_burn_in + 1):n_days

# Seasonality
x_eta = 0.00:0.002:0.5
length(x_eta)

# Waning
rho_step = 0.00002
x_rho = rho_step:rho_step:0.005
length(x_rho)

x_vals = vec([(x1, x2) for x1 in x_eta, x2 in x_rho])

ix_jobs = get_jobs(arg_ix, n_array, length(x_vals))
x_vals_job = x_vals[ix_jobs]



y_period = zeros(length(x_vals_job), 3)
y_inf_maxima = zeros(length(x_vals_job), 2)

time_start = Base.time()

Threads.@threads for i in eachindex(x_vals_job)
    println("Running job $(ix_jobs[i])")

    model_params = make_model_parameters(
        k = k, beta = beta, gamma = gamma, C = C, rho = x_vals_job[i][2],
        b = b, h = h, c_jump_dist = c_jump_dist; boosting = "none", eta = x_vals_job[i][1]
    )

    ode_solution = ode_solve(model_params, n_days, n_inf_0, ode_sparsity, saveat = Δt)

    y = ode_solution(t)[1:(model_params.S * 2), :]'

    y_daily = ode_solution(t_daily)[(model_params.S + 1):(model_params.S * 2), :]'
    y_inf_maxima[i, 1] = minimum(sum(y_daily, dims = 2))
    y_inf_maxima[i, 2] = maximum(sum(y_daily, dims = 2))

    y_zeroed = copy(y)

    for j in axes(y, 1)
        y_zeroed[j, :] .= y[j, :] .- y[1, :]
    end
    sum_abs2 = vec(sum(abs2.(y_zeroed), dims = 2))

    t_repeats = t[findall(sum_abs2 .< ϵ)]
    n_repeats = length(t_repeats)
    t_repeat_groups = zeros(Int, length(t_repeats))
    t_repeat_groups[1] = 1

    for i in 2:n_repeats
        diff = t_repeats[i] - t_repeats[i - 1]
        if diff > 10.0
            t_repeat_groups[i] = t_repeat_groups[i - 1] + 1
        else
            t_repeat_groups[i] = t_repeat_groups[i - 1]
        end
    end

    t_repeats_grouped = zeros(maximum(t_repeat_groups))

    for i in eachindex(t_repeats_grouped)
        t_repeats_grouped[i] = mean(t_repeats[t_repeat_groups .== i])
    end

    period_samples = (t_repeats_grouped .- t_repeats_grouped[1]) ./ (0:(length(t_repeats_grouped) - 1))
    y_period[i, 1] = mean(period_samples[2:end])
    y_period[i, 2] = std(period_samples[2:end])
    y_period[i, 3] = length(period_samples[2:end])
end

time_elapsed = Base.time() - time_start

println("Completed $(length(ix_jobs)) jobs in $(round(time_elapsed, digits = 2)), ($(round(time_elapsed/length(ix_jobs), digits = 2)) seconds/job)")

t_seq = collect(t_daily)
x_vals_job = stack(x_vals_job)
jldsave("data_dist/period_grid/$(arg_ix).jld2"; x_vals_job, t_seq, y_period, y_inf_maxima)

println("Outputs saved.")