

include("dependencies.jl")

using ForwardDiff
using DataFrames

c_max = 8
k = 4
#gamma = 0.5
beta = 1.5
b = 3
m = 4
c_jump_dist = Normal(7, 0.5)


x_gamma = collect(0.1:0.005:1.0)
x_lambda = collect(0.01:0.005:0.4)

x_df = expand_grid(;x_gamma, x_lambda)

y_eigvals = Vector{Vector{Complex{Double64}}}(undef, nrow(x_df))


Threads.@threads for (i, r) in collect(enumerate(eachrow(x_df)))
    println(i)

    model_params = make_model_parameters(
        c_max = c_max, k = k, beta = beta, gamma = r.x_gamma, lambda = r.x_lambda,
        b = b, m = m, c_jump_dist = c_jump_dist
    )
    
    sus_steady, inf_steady = get_steady_state(model_params)
    

    u0 = zeros(Double64, n_compartments * model_params.N)
    u0[ode_ix(c_sus, 1:model_params.N, model_params.N)] = sus_steady
    u0[ode_ix(c_inf, 1:model_params.N, model_params.N)] = inf_steady
    
    
    du0 = copy(u0)
    J = ForwardDiff.jacobian((du, u) -> ode_step_no_count!(du, u, model_params, 0.0), du0, u0)
    y_eigvals[i] = eigvals(J)
end

mat_real_eigvals = reduce(vcat,transpose.(real(y_eigvals)))
vec_real_eigvals = maximum(mat_real_eigvals, dims = 2)

x_df[!, "max_eigval"] .= vec_real_eigvals

using CSV

CSV.write("data/fractal.csv", x_df)


