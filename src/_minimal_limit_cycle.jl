
include("dependencies.jl")

k = 6

R = 1.5
gi = 0.5

model_params = make_model_parameters(
    c_max = 1, k = k,

    beta = R / gi, gamma = 1 / gi, lambda = 0.02,

    b = 1 / k + 1e-5, m = Inf,

    c_jump_dist = Dirac(2)
)



heatmap(model_params.M)
plot(model_params.c_levels, model_params.p_acq, seriestype = :scatter, legend = false)

ode_sparsity = ode_get_sparsity(model_params)

n_inf_0 = 0.001
n_days = 1000
saveat = 1

ode_solution = ode_solve(
    model_params, n_days, n_inf_0, ode_sparsity;
    saveat = saveat, dt = 0.01, datatype = Double64
);

plot(log.(ode_solution)', legend = false)

susceptible = zeros(BigFloat, n_days ÷ saveat, model_params.N)
infected = zeros(BigFloat, n_days ÷ saveat, model_params.N)
counts = zeros(BigFloat, n_days ÷ saveat, model_params.N)

for k in 1:model_params.N
    for d in 1:(n_days ÷ saveat)
        susceptible[d, k] = ode_solution(d * saveat)[ode_ix(c_sus, k, model_params.N)]
        infected[d, k] = ode_solution(d * saveat)[ode_ix(c_inf, k, model_params.N)]
        counts[d, k] = ode_solution(d * saveat)[ode_ix(c_count, k, model_params.N)]
    end
end


plot(infected)

plot(log.(infected))
plot(susceptible)
plot(log.(susceptible))
plot((sum(infected, dims = 2)))
plot(log.(sum(infected, dims = 2)))

heatmap(min.(susceptible[1:end,:]',0.1))


