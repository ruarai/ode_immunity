
include("src/dependencies.jl")

model_params = make_model_parameters(
    c_max = 8, k = 4,

    beta = 1.5, gamma = 0.5, lambda = 0.1,

    b = 3, m = 4,

    c_jump_dist = Normal(7, 0.5)
)


ode_sparsity = ode_get_sparsity(model_params)

n_inf_0 = 0.01
n_days = 8000

ode_solution = @time ode_solve(model_params, n_days, n_inf_0, ode_sparsity)

infected = zeros(n_days, model_params.N)

for k in 1:model_params.N
    for d in 1:n_days
        infected[d, k] = sum(ode_solution(d)[ode_ix(c_inf, k, model_params.N)])
    end
end


plot((sum(infected, dims = 2)))



u0_steady = get_steady_state(model_params)





ode_step_fn = ODEFunction(ode_step_no_count!; jac_prototype = float.(ode_sparsity))

u0 = zeros(Double64, n_compartments * model_params.N)
u0[ode_ix(c_sus, 1:model_params.N, model_params.N)] = u0_steady[1:model_params.N]
u0[ode_ix(c_inf, 1:model_params.N, model_params.N)] = u0_steady[(1:model_params.N) .+ model_params.N]

du = zeros(Double64, n_compartments * model_params.N)
ode_step_no_count!(du, u0, model_params, 0.0)
du

n_days = 1000

tspan = (0.0, 1.0 * n_days)
prob = ODEProblem{true, SciMLBase.FullSpecialize}(ode_step_fn, u0, tspan, model_params)

sol_steady = DifferentialEquations.solve(prob, Euler(), dt = 0.01, saveat = 1)


plot(sol_steady, legend = false)



omega_inv = 1 ./ (1 .- model_params.p_acq)

N = model_params.N

u0 = SVector{N}(convert.(Double64, [1.0 / N for x in 1:N]))

fn_solve(inf_vec, p) = steady_state_and_valid(
    inf_vec, 
    model_params.B, model_params.M, 
    omega_inv, 
    model_params.k, model_params.lambda, model_params.beta, model_params.gamma
)

probN = NonlinearProblem(fn_solve, u0)
inf_vec = solve(
    probN, NewtonRaphson();
    abstol = 1e-50, maxiters = 200,
    show_trace = Val(true), trace_level = NonlinearSolve.TraceAll(10)
)


sus_vec = (model_params.gamma / (model_params.beta * sum(inf_vec))) * (inf_vec .* omega_inv)

I_sum = sum(inf_vec)

a = steady_state(inf_vec, I_sum, model_params.B, model_params.M,
                 omega_inv, model_params.k, model_params.lambda, model_params.beta)
b = 1 - population_sum(inf_vec, I_sum, omega_inv, model_params.beta, model_params.gamma)
c = soft_inf_positive(inf_vec)