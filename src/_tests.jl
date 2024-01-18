
include("dependencies.jl")

model_params = make_model_parameters(
    c_max = 8, k = 4,

    beta = 1.5, gamma = 0.5, lambda = 0.4,

    b = 3, m = 4,

    c_jump_dist = Normal(7, 0.5)
)


ode_sparsity = ode_get_sparsity(model_params)

n_inf_0 = 0.01
n_days = 2000

ode_solution = @time ode_solve(model_params, n_days, n_inf_0, ode_sparsity)

infected = zeros(n_days, model_params.N)

for k in 1:model_params.N
    for d in 1:n_days
        infected[d, k] = sum(ode_solution(d)[ode_ix(c_inf, k, model_params.N)])
    end
end


plot((sum(infected, dims = 2)))



sus_steady, inf_steady = get_steady_state(model_params)





ode_step_fn = ODEFunction(ode_step_no_count!; jac_prototype = float.(ode_sparsity))

u0 = zeros(Double64, n_compartments * model_params.N)
u0[ode_ix(c_sus, 1:model_params.N, model_params.N)] = sus_steady
u0[ode_ix(c_inf, 1:model_params.N, model_params.N)] = inf_steady

du = zeros(Double64, n_compartments * model_params.N)
ode_step_no_count!(du, u0, model_params, 0.0)
maximum(du)

n_days = 1000

tspan = (0.0, 1.0 * n_days)
prob = ODEProblem{true, SciMLBase.FullSpecialize}(ode_step_fn, u0, tspan, model_params)

sol_steady = DifferentialEquations.solve(prob, Euler(), dt = 0.01, saveat = 1)


plot(sol_steady, legend = false)


using ForwardDiff
    
du0 = copy(u0)
J = ForwardDiff.jacobian((du, u) -> ode_step_no_count!(du, u, model_params, 0.0), du0, u0)
J_eigvals = eigvals(J, maxiter = 1000, reltol = 1e-30)

plot(J_eigvals, seriestype = :scatter, xlim = (-0.1, 0.1))


println(real(J_eigvals)[end])


