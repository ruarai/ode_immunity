
using JLD2


jld_files = readdir("data_dist/period_grid", join = true)

t_seq = load(jld_files[1], "t_seq")


n_files = length(jld_files)

x_vals = Vector{Matrix{Float64}}(undef, n_files)
y_inf_maxima = Vector{Matrix{Float64}}(undef, n_files)
y_period = Vector{Matrix{Float64}}(undef, n_files)


for i in 1:n_files
    print("$i ")
    f = jld_files[i]

    x_vals[i] = load(f, "x_vals_job")
    y_inf_maxima[i] = load(f, "y_inf_maxima")
    y_period[i] = load(f, "y_period")
end


x_vals = reduce(hcat, x_vals)
y_inf_maxima = reduce(vcat, y_inf_maxima)
y_period = reduce(vcat, y_period)

jldsave("data/paper/period_over_grid.jld2"; x_vals, y_inf_maxima, y_period)
