include("dependencies.jl")

using JLD2

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
    b = b, h = h, c_jump_dist = c_jump_dist; boosting = "none"
)
ode_sparsity = ode_get_sparsity(model_params_0)

n_inf_0 = 0.01
n_days = 32000



x_eta = 0:0.01:1


y_fixed_I = zeros(length(x_eta))

y_I_sol = zeros(length(x_eta), n_days)
y_count_sol = zeros(length(x_eta), n_days)

Threads.@threads for i in eachindex(x_eta)

    model_params = make_model_parameters(
        k = k, beta = beta, gamma = gamma, C = C, rho = rho,
        b = b, h = h, c_jump_dist = c_jump_dist; eta = x_eta[i], boosting = "none"
    )

    # Calculate the fixed point/steady state solution
    # sus_vec, inf_vec = get_steady_state(model_params)
    # y_fixed_I[i] = sum(inf_vec)


    # Calculate a (not-necessarily-stable) solution of the ODE over time
    ode_solution = ode_solve(model_params, n_days, n_inf_0, ode_sparsity)

    for d in 1:n_days
        y_I_sol[i, d] = sum(ode_solution(d)[ode_ix(c_inf, 1:model_params.S, model_params.S)])
        y_count_sol[i, d] = sum(ode_solution(d)[ode_ix(c_count, 1:model_params.S, model_params.S)])
    end

    println(i)
end

plot(maximum(y_I_sol[:, 30000:end], dims = 2))
plot!(x_eta, minimum(y_I_sol[:, 30000:end], dims = 2), yaxis = :log10)


plot(y_I_sol[35, 25000:end], ylim = (0, 0.2))

plot(y_I_sol[6, 25000:end])
vline!([i * 365 for i in 0:20])
plot(y_I_sol[3, 25000:end])

heatmap(log10.(y_I_sol[:,1:8000]))
heatmap(log10.(y_I_sol[:,30000:end]))

plot(y_I_sol[:,end], legend = false)


y_inc_sol = diff(y_count_sol, dims = 2)

burn_in_days = 30000


period = zeros(Int, length(x_eta))
attack_rate = zeros(length(x_eta))

@showprogress Threads.@threads for i in eachindex(x_eta)
    incidence = y_inc_sol[i, burn_in_days:end]
    period[i] = get_period(incidence)
    attack_rate[i] = get_periodic_attack_rate(incidence, period[i])
end


plot(period)
hline!([365])
plot(attack_rate)
vline!([75])

plot(y_inc_sol[75, burn_in_days:end])
vline!([i*365 for i in 0:3])

get_period(y_inc_sol[75, burn_in_days:end])



y = y_inc_sol[25, burn_in_days:end]

using FFTW

plot(y)

fft_y = fft(y)
fft_freq_y = fftfreq(length(y), 1) .* length(y)


plot(fft_freq_y, abs.(fft_y))

n = length(y)
ns = isodd(n) ? (n-1)÷2 : n÷2-1 
AYf = abs.(fft_y)[1:ns]
fft_freq_y = fft_freq_y[1:ns]

plot(length(y) ./ fft_freq_y, AYf)

plot(1 ./ fft_freq_y, imag.(fft_y))


ix = partialsortperm(AYf, 1:2, rev=true)
f1, f2 = fft_freq_y[ix]
Δf = abs(f2 - f1)   
T₀ = 2/Δf    