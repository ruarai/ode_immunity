include("dependencies.jl")

function decay_step!(du, u, p, t)
    du[1] = -p.r * u[1]
end

u0 = [10.0^8]
n_days = 100

function decay_strata!(du, u, p, t)
    decay_rate = p.lambda * p.k

    du[1] = u[2] * decay_rate
    for i in 2:p.k
        du[i] = u[i + 1] * decay_rate - u[i] * decay_rate
    end
    du[k + 1] = -u[k + 1] * decay_rate
end


k = 32
C = 8
lambda = 0.01

c_levels = 10.0 .^ (C * (0:k) ./ k)

# lambda = r / (C * log(10))

# lambda = r * (k * n - 1) / (C * log(10))

strata_params = (k = k, lambda = lambda)


u0_strata = zeros(k + 1)
u0_strata[k+1] = 1.0

du0 = copy(u0_strata)
ode_sparsity = Symbolics.jacobian_sparsity((du, u) -> decay_strata!(du, u, strata_params, 0.0), du0, u0_strata)
strata_step_fn = ODEFunction(decay_strata!; jac_prototype = float.(ode_sparsity))
prob_strata = ODEProblem(decay_strata!, u0_strata, (0.0, n_days), strata_params)

sol_strata = DifferentialEquations.solve(prob_strata, Euler(), dt = 0.01, saveat = 0.5)

y_strata = stack([sol_strata(i) for i in 1:n_days])
heatmap(1:n_days, log10.(c_levels), min.(0.1, y_strata))


avg = [sum(y_strata[:, t] .* c_levels) for t in 1:n_days]


r = k * lambda * (1 - 10^(-C/k))



prob_static = ODEProblem(decay_step!, u0, (0.0, n_days), (r = r, ))
sol_static = DifferentialEquations.solve(prob_static, Euler(), dt = 0.01, saveat = 0.5)

heatmap(1:n_days, log10.(c_levels), min.(0.2, y_strata))
plot!(log10.(avg))
plot!(log10.(sol_static(1:n_days))', lc = "red")
