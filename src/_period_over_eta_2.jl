include("dependencies.jl")

k = 32
C = 8.0
R = 1.5
gamma = 0.25
beta = R * gamma
rho = 0.003
b = 2^3
h = 8
c_jump_dist = Normal(2^6, 2^5)


model_params_0 = make_model_parameters(
    k = k, beta = beta, gamma = gamma, C = C, rho = rho,
    b = b, h = h, c_jump_dist = c_jump_dist; boosting = "none", eta = 0.0
)

ode_sparsity = ode_get_sparsity(model_params_0)


n_inf_0 = 0.0001
n_days = 60000
n_days_burn_in = 30000

ϵ = 1e-6
Δt = 0.25
t = n_days_burn_in:Δt:n_days
t_daily = (n_days_burn_in + 1):n_days

x_eta = 0:0.001:0.5
y_period = zeros(length(x_eta), 3)
y_inf = zeros(length(x_eta), n_days - n_days_burn_in)


@showprogress Threads.@threads for i_eta in eachindex(x_eta)

    model_params = make_model_parameters(
        k = k, beta = beta, gamma = gamma, C = C, rho = rho,
        b = b, h = h, c_jump_dist = c_jump_dist; boosting = "none", eta = x_eta[i_eta]
    )

    ode_solution = ode_solve(model_params, n_days, n_inf_0, ode_sparsity, saveat = Δt)

    y = ode_solution(t)[1:(model_params.S * 2), :]'

    y_daily = ode_solution(t_daily)[(model_params.S + 1):(model_params.S * 2), :]'
    y_inf[i_eta, :] = vec(sum(y_daily, dims = 2))

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
    y_period[i_eta, 1] = mean(period_samples[2:end])
    y_period[i_eta, 2] = std(period_samples[2:end])
    y_period[i_eta, 3] = length(period_samples[2:end])
end

scatter(x_eta, y_period[:, 2])
scatter(x_eta, y_period[:, 3] .== 1)
good_est = (y_period[:, 2] .< 1.0) .& (y_period[:, 3] .> 1.0)

scatter(x_eta[good_est], y_period[good_est, 1])


scatter(x_eta[good_est], y_period[good_est, 1] ./ 365)

heatmap( log.(y_inf))
heatmap((n_days - 10000):n_days, x_eta, log.(y_inf[:,(n_days - 10000):n_days]))


i = 48
scatter(x_eta[good_est], y_period[good_est, 1], ylim = (0, 2000))
vline!([x_eta[i]])
plot(log.(y_inf[i,1:2000]))
plot!(log.(y_inf[i, round(Int, y_period[i, 1]):2000]))
vline!([j * y_period[i, 1] for j in 1:2])