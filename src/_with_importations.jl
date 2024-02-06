
include("dependencies.jl")

c_max = 8
k = 2
R = 1.5
gamma = 0.25
beta = R * gamma
lambda = 0.01
b = 2
m = 5
c_jump_dist = Normal(4.5, 0.5)

model_params = make_model_parameters(
    c_max = c_max, k = k, beta = beta, gamma = gamma, lambda = lambda,
    b = b, m = m, c_jump_dist = c_jump_dist
)


function ode_step_import!(du, u, model_params, t)
    N = model_params.N

    for i in 1:(N * n_compartments)
        du[i] = 0
    end

    flow_decay = @views u[ode_ix(1, 2:N, N)] .* model_params.lambda .* model_params.k

    du[ode_ix(1, 2:N, N)] = -flow_decay
    du[ode_ix(1, 1:(N - 1), N)] .+= flow_decay

    flow_sus_to_inf = @views u[ode_ix(c_sus, 1:N, N)] .*
        sum(u[ode_ix(c_inf, 1:N, N)]) .*
        model_params.beta .* 
        (1 .- model_params.p_acq) .+

        u[ode_ix(c_sus, 1:N, N)] *
        0.000001

    flow_inf_to_sus = @views u[ode_ix(c_inf, 1:N, N)] .* model_params.gamma


    du[ode_ix(c_sus, 1:N, N)] .+= model_params.M * flow_inf_to_sus - flow_sus_to_inf
    du[ode_ix(c_inf, 1:N, N)] .+= flow_sus_to_inf - flow_inf_to_sus

    du[ode_ix(c_count, 1:N, N)] .+= flow_sus_to_inf 
end


ode_step_fn = ODEFunction(ode_step_import!)

u0 = zeros(Float64, n_compartments * model_params.N)
u0[ode_ix(c_sus, 1, model_params.N)] = 1.0 - n_inf_0
u0[ode_ix(c_inf, 1, model_params.N)] = n_inf_0

n_days = 10000
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

plot(infected)

plot(log.(infected))
plot(susceptible)
plot(log.(susceptible))
plot((sum(infected, dims = 2)))
plot(log.(sum(infected, dims = 2)))

heatmap(min.(susceptible',0.1))