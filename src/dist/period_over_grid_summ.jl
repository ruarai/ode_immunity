
using JLD2


jld_files = readdir("data_dist/period_grid", join = true)

n_files = length(jld_files)

x_vals = Vector{Matrix{Float64}}(undef, n_files)
y_inf_summary = Vector{Matrix{Float64}}(undef, n_files)
y_period = Vector{Matrix{Float64}}(undef, n_files)
y_attack_rate = Vector{Vector{Float64}}(undef, n_files)

complete_tasks = Vector{Int}(undef, 0)


for i in 1:n_files
    print("$i ")
    data = load(jld_files[i])

    x_vals[i] = data["x_vals_job"]
    y_inf_summary[i] = data["y_inf_summary"]
    y_period[i] = data["y_period"]
    y_attack_rate[i] = data["y_attack_rate"]

    append!(complete_tasks, parse(Int, match(r"\d+", jld_files[i]).match))
end

incomplete_tasks = setdiff(1:8350, complete_tasks)'
incomplete_task_string = join(incomplete_tasks, ',')


x_vals = reduce(hcat, x_vals)
y_inf_summary = reduce(vcat, y_inf_summary)
y_period = reduce(vcat, y_period)

y_attack_rate = reduce(vcat, y_attack_rate)

jldsave("data/paper/period_over_grid.jld2"; x_vals, y_inf_summary, y_period, y_attack_rate)
