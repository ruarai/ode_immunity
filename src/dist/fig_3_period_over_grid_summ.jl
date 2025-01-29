using JLD2


jld_files = readdir("data_dist/period_grid", join = true)

n_files = length(jld_files)

x_vals = Vector{Matrix{Float64}}(undef, 0)
y_inf_summary = Vector{Matrix{Float64}}(undef, 0)
y_period = Vector{Matrix{Float64}}(undef, 0)
y_seasonality = Vector{Matrix{Float64}}(undef, 0)

y_peaks = Vector{Tuple{Float64, Float64, Float64, Float64}}(undef, 0)

complete_tasks = Vector{Int}(undef, 0)



for i in 1:n_files
    print("$i ")
    data = load(jld_files[i])

    if !haskey(data, "y_seasonality")
        continue
    end

    push!(x_vals, data["x_vals_job"])
    push!(y_inf_summary, data["y_inf_summary"])
    push!(y_period, data["y_period"])
    push!(y_seasonality, data["y_seasonality"])

    # y_peaks_vec = data["y_peaks"]
    # for j in eachindex(y_peaks_vec)
    #     for k in eachindex(y_peaks_vec[j])
    #         push!(y_peaks, (x_vals[i][1, j], x_vals[i][2, j], y_peaks_vec[j][k][1], y_peaks_vec[j][k][2]))
    #     end
    # end

    append!(complete_tasks, parse(Int, match(r"\d+", jld_files[i]).match))
end

incomplete_tasks = setdiff(1:8350, complete_tasks)'
incomplete_task_string = join(incomplete_tasks, ',')


x_vals = reduce(hcat, x_vals)
y_inf_summary = reduce(vcat, y_inf_summary)
y_period = reduce(vcat, y_period)
y_seasonality = reduce(vcat, y_seasonality)

y_peaks = Matrix(reinterpret(reshape, Float64, y_peaks))


jldsave("data/paper/period_over_grid.jld2"; x_vals, y_inf_summary, y_period, y_seasonality, y_peaks)
