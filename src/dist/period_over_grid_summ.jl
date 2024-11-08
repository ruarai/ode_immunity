
using JLD2


jld_files = readdir("data_dist/period_grid", join = true)

t_seq = load(jld_files[1], "t_seq")


n_files = length(jld_files)

x_vals = Vector{Matrix{Float64}}(undef, n_files)
y_inf_max = Vector{Vector{Float64}}(undef, n_files)
y_inf_min = Vector{Vector{Float64}}(undef, n_files)

for i in 1:n_files
    f = jld_files[i]

    x_vals[i] = load(f, "x_vals_job")

    y_inf = load(f, "y_inf")
    y_inf_max[i] = vec(maximum(y_inf, dims = 2))
    y_inf_min[i] = vec(maximum(y_inf, dims = 2))
end


x_vals = reduce(hcat, x_vals)
y_inf_max = reduce(vcat, y_inf_max)
y_inf_min = reduce(vcat, y_inf_min)

jldsave("data/paper/period_over_grid.jld2"; x_vals, y_inf_max, y_inf_min)