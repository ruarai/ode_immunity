

include("dependencies.jl")

# results_peak_mag = Array{Array{Float64}}(undef, n_params)
# results_troughs_mag = Array{Array{Float64}}(undef, n_params)

# results_peak_trough_diff = Array{Array{Float64}}(undef, n_params)

# results_peak_mag_diff_final = Array{Float64}(undef, n_params)
# results_peak_trough_diff_final = Array{Float64}(undef, n_params)

# for i in 1:n_params
#     res = results[i, :]
#     res_diff = diff(res)
#     res_peaks = res_diff[1:(end - 1)] .> 0 .&& res_diff[2:end] .<= 0
#     results_peak_mag[i] = results[i, 2:(end - 1)][res_peaks]

#     res_troughs = res_diff[1:(end - 1)] .< 0 .&& res_diff[2:end] .>= 0
#     results_troughs_mag[i] = results[i, 2:(end - 1)][res_troughs]

#     min_common_len = min(length(results_peak_mag[i]), length(results_troughs_mag[i]))

#     results_peak_trough_diff[i] = zeros(min_common_len)

#     for j in 1:min_common_len
#         results_peak_trough_diff[i][j] = results_peak_mag[i][j] - results_troughs_mag[i][j]
#     end

#     if length(results_peak_mag[i]) > 1
#         results_peak_mag_diff_final[i] = diff(results_peak_mag[i])[end]
#         results_peak_trough_diff_final[i] = results_peak_trough_diff[i][end]
#     end
# end

using DataFrames
using JLD2
results = load("data/paper/delay_de.jld2")

plot(results["results_peak_mag"], legend = false)
plot!(results["results_troughs_mag"], legend = false)

plot(results["results_peak_trough_diff"], legend = false)
plot(results_peak_trough_diff_final .> 0.01)


heatmap(
    R0s, lambdas,
    log.(reshape(results["results_peak_trough_diff_final"], length(R0s), length(lambdas)) .+ 0.01)' .> log(0.01) + 1e-2
)
heatmap(
    R0s, lambdas,
    log.(reshape(results["results_peak_trough_diff_final"], length(R0s), length(lambdas)) .+ 0.01)'
)


taus = 1.0 ./ lambdas

taus * 0.5

2 .* taus .* 0.5 .+ 3


R0 <= 2 τ * γ + 3