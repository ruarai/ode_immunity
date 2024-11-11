
using JLD2


jld_files = readdir("data_dist/period_grid", join = true)

n_files = length(jld_files)

x_vals = Vector{Matrix{Float64}}(undef, n_files)
y_inf_maxima = Vector{Matrix{Float64}}(undef, n_files)
y_period = Vector{Matrix{Float64}}(undef, n_files)
y_attack_rate = Vector{Vector{Float64}}(undef, n_files)


for i in 1:n_files
    print("$i ")
    data = load(jld_files[i])

    x_vals[i] = data["x_vals_job"]
    y_inf_maxima[i] = data["y_inf_maxima"]
    y_period[i] = data["y_period"]

    if haskey(data, "y_attack_rate")
        y_attack_rate[i] = data["y_attack_rate"]
    else
        y_attack_rate[i] = zeros(size(data["x_vals_job"])[2]) .* NaN
    end
end


x_vals = reduce(hcat, x_vals)
y_inf_maxima = reduce(vcat, y_inf_maxima)
y_period = reduce(vcat, y_period)

y_attack_rate = reduce(vcat, y_attack_rate)

jldsave("data/paper/period_over_grid.jld2"; x_vals, y_inf_maxima, y_period, y_attack_rate)
