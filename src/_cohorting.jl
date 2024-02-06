
include("dependencies.jl")

c_max = 1
k = 100
R = 1.8
gamma = 0.25
beta = R * gamma
lambda = 0.005
b = 1 / k + 1e-10
m = Inf
c_jump_dist = Normal(1, 1e-10)

model_params = make_model_parameters(
    c_max = c_max, k = k, beta = beta, gamma = gamma, lambda = lambda,
    b = b, m = m, c_jump_dist = c_jump_dist,

    boosting = false
)

plot(model_params.c_levels, model_params.p_acq, seriestype = :scatter)
heatmap(model_params.M)

function ode_step_beta!(du, u, model_params, t)
    N = model_params.N

    for i in 1:(N * n_compartments)
        du[i] = 0
    end

    flow_decay = @views u[ode_ix(1, 2:N, N)] .* model_params.lambda .* model_params.k

    du[ode_ix(1, 2:N, N)] = -flow_decay
    du[ode_ix(1, 1:(N - 1), N)] .+= flow_decay

    t_mult = 1.0

    if t > 60
        t_mult = 1.0
    elseif t > 25
        t_mult = 0.3
    end

    beta_t = model_params.beta * t_mult

    flow_sus_to_inf = @views u[ode_ix(c_sus, 1:N, N)] .*
        sum(u[ode_ix(c_inf, 1:N, N)]) .*
        beta_t .* 
        (1 .- model_params.p_acq)

    flow_inf_to_sus = @views u[ode_ix(c_inf, 1:N, N)] .* model_params.gamma


    du[ode_ix(c_sus, 1:N, N)] .+= model_params.M * flow_inf_to_sus - flow_sus_to_inf
    du[ode_ix(c_inf, 1:N, N)] .+= flow_sus_to_inf - flow_inf_to_sus

    du[ode_ix(c_count, 1:N, N)] .+= flow_sus_to_inf 
end


ode_sparsity = ode_get_sparsity(model_params)

ode_step_fn = ODEFunction(ode_step_beta!; jac_prototype = float.(ode_sparsity))

u0 = zeros(Float64, n_compartments * model_params.N)
u0[ode_ix(c_sus, 1, model_params.N)] = 1.0 - n_inf_0
u0[ode_ix(c_inf, 1, model_params.N)] = n_inf_0

n_days = 4000
tspan = (0.0, 1.0 * n_days)
prob = ODEProblem{true, SciMLBase.FullSpecialize}(ode_step_fn, u0, tspan, model_params)

sol = DifferentialEquations.solve(prob, Euler(), dt = 0.01, saveat = 1);


susceptible = zeros(BigFloat, n_days, model_params.N)
infected = zeros(BigFloat, n_days, model_params.N)

for k in 1:model_params.N
    for d in 1:(n_days)
        susceptible[d, k] = sol(d)[ode_ix(c_sus, k, model_params.N)]
        infected[d, k] = sol(d)[ode_ix(c_inf, k, model_params.N)]
    end
end

plot(infected, legend = false)

plot(log.(infected), legend = false)
plot(susceptible, legend = false)
plot(log.(susceptible))
plot((sum(infected, dims = 2)))

heatmap(min.(susceptible',0.01))

plot(log.(sum(infected, dims = 2)))
plot(log.(sum(susceptible, dims = 2)))