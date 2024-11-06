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



n_inf_0 = 0.0001
n_days = 120000
n_days_burn_in = 40000

Δt = 0.1
t = n_days_burn_in:Δt:n_days

model_params = make_model_parameters(
    k = k, beta = beta, gamma = gamma, C = C, rho = rho,
    b = b, h = h, c_jump_dist = c_jump_dist; boosting = "none", eta = 0.1
)

ode_sparsity = ode_get_sparsity(model_params)
ode_solution = ode_solve(model_params, n_days, n_inf_0, ode_sparsity, saveat = 0.1)


y = ode_solution(t)[1:(model_params.S * 2), :]'
y_zeroed = copy(y)

for j in axes(y, 1)
    y_zeroed[j, :] .= y[j, :] .- y[1, :]
end

plot(y_zeroed[1:1:6000,:], lc = "black", legend = false)

sum_abs2 = vec(sum(abs2.(y_zeroed), dims = 2))


ϵ = 1e-6
plot(sum_abs2[2:end], yaxis = :log10)
plot(sum_abs2[1:end])

plot((sum_abs2 .< ϵ)[1:end])

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

plot(t_repeats_grouped)

period_samples = (t_repeats_grouped .- t_repeats_grouped[1]) ./ (0:(length(t_repeats_grouped) - 1))


histogram(period_samples)
# If std is high, our estimate is bad
std(period_samples[2:end])
period = median(period_samples[2:end])

plot(vec(sum(y[:, (model_params.S + 1):(model_params.S * 2)], dims = 2))[1:40000], legend = false)


(period) / 365
period_ix = round(Int, period / Δt)

plot(y_zeroed[1:10:60000,:], lc = "black", legend = false)
plot!(y_zeroed[period_ix:10:(period_ix + 60000),:], lc = "red", alpha = 0.5, legend = false)




vline!([i * period / Δt for i in 1:10])

