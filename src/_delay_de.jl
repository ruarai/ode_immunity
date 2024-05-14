

include("dependencies.jl")

using DataFrames
using JLD2


function sim_dde(; R0, γ, τ, dt, n_days, n_inf_0)
    t_steps = 0:dt:n_days
    n_t_steps = length(t_steps)

    β = R0 * γ

    S_t = zeros(n_t_steps)
    I_t = zeros(n_t_steps)
    R_t = zeros(n_t_steps)
    dI_t = zeros(n_t_steps)

    S_t[1] = 1 - n_inf_0
    I_t[1] = n_inf_0

    P_t = 1 .- cdf.(τ, t_steps)
    ix_P_t_neglible = findfirst(P_t .< 1e-8)
    if isnothing(ix_P_t_neglible)
        ix_P_t_neglible = length(P_t)
    end
    
    for t in 1:(n_t_steps - 1)
        I = I_t[t]
    
        recovered = 0
        if t > 1
            for x in 1:min(t - 1, ix_P_t_neglible)
                recovered += @views I_t[t - x] * P_t[x + 1] * dt
            end
        end

        R = γ * recovered

    
        S = 1 - I - R
    
        dI = -γ * I + β * I * S
    
        S_t[t + 1] = S - dI * dt
        I_t[t + 1] = I + dI * dt
        R_t[t + 1] = 1 -  S_t[t + 1] - I_t[t + 1]
        dI_t[t + 1] = dI 

    end


    return(0:dt:n_days, S_t, I_t, R_t, dI_t)

    # return (0:n_days,
    #         S_t[1:round(Int, 1 / dt):end],
    #         I_t[1:round(Int, 1 / dt):end],
    #         R_t[1:round(Int, 1 / dt):end],
    #         dI_t[1:round(Int, 1 / dt):end])
end


n_days = 365*10


println("Initial run.")

c_jump_dist = Normal(0.8, 0.05)
lambda = 0.002

τ = c_jump_dist / lambda

t_days, S_t, I_t, R_t, dI_t = sim_dde(;
    R0 = 1.5,
    γ = 0.25,
    τ = τ,
    dt = 0.01,
    n_days, 
    n_inf_0 = 0.001
)

plot(t_days, I_t)
plot(t_days, S_t)
plot(t_days, R_t)
plot(t_days, dI_t)

plot(t_days, S_t .+ I_t .+ R_t)


lambdas = 0.001:0.001:0.05
R0s = 1.0:0.02:3.5

param_sets = expand_grid(R0 = R0s, lambda = lambdas)

n_params = size(param_sets, 1)
results_peak_mag = Array{Array{Float64}}(undef, n_params)
results_troughs_mag = Array{Array{Float64}}(undef, n_params)

results_peak_trough_diff = Array{Array{Float64}}(undef, n_params)

results_peak_mag_diff_final = Array{Float64}(undef, n_params)
results_peak_trough_diff_final = Array{Float64}(undef, n_params)


println("Running across $n_params sets of parameters.")

using ProgressMeter
p = Progress(n_params)

Threads.@threads for i in axes(param_sets, 1)
    row = param_sets[i, :]

    τ = c_jump_dist / row.lambda
    t_days, S_t, I_t, R_t, dI_t = sim_dde(;
        R0 = row.R0, 
        γ = 0.5,
        τ = τ,
        dt = 0.01,
        n_days, 
        n_inf_0 = 0.01
    )

    res_diff = diff(I_t)
    res_peaks = res_diff[1:(end - 1)] .> 0 .&& res_diff[2:end] .<= 0
    results_peak_mag[i] = I_t[2:(end - 1)][res_peaks]

    res_troughs = res_diff[1:(end - 1)] .< 0 .&& res_diff[2:end] .>= 0
    results_troughs_mag[i] = I_t[2:(end - 1)][res_troughs]

    min_common_len = min(length(results_peak_mag[i]), length(results_troughs_mag[i]))

    results_peak_trough_diff[i] = zeros(min_common_len)

    for j in 1:min_common_len
        results_peak_trough_diff[i][j] = results_peak_mag[i][j] - results_troughs_mag[i][j]
    end

    if length(results_peak_mag[i]) > 1
        results_peak_mag_diff_final[i] = diff(results_peak_mag[i])[end]
        results_peak_trough_diff_final[i] = results_peak_trough_diff[i][end]
    end

    next!(p)
end
finish!(p)

jldsave("data/paper/delay_de.jld2"; 
    param_sets,
    results_peak_mag,
    results_troughs_mag,
    results_peak_trough_diff,
    results_peak_mag_diff_final,
    results_peak_trough_diff_final
)

