include("dependencies.jl")


using MultiFloats
using Decimals

function ode_step_sirs!(du, u, sirs_params, t)
   
    S = u[1]
    I = u[2]
    R = u[3]

    β, γ, δ = sirs_params 

    du[1] = - β * S * I + δ * R
    du[2] = β * S * I - γ * I
    du[3] = γ * I - δ * R
end

n_inf_0 = 0.001

u0 = zeros(Double64, 3)
u0[1] = 1 - n_inf_0
u0[2] = n_inf_0

sirs_params = (2.0, 1.999, 0.001)
sirs_params[1] / sirs_params[2]


ode_step_fn = ODEFunction(ode_step_sirs!)


n_days = 100000
tspan = (0.0, 1.0 * n_days)
prob = ODEProblem{true, SciMLBase.FullSpecialize}(ode_step_fn, u0, tspan, sirs_params)

sol = DifferentialEquations.solve(
    prob, Euler(),
    dt = 0.01, saveat = 10.0
);

plot(sol[:,1:2:end])
plot(log.(sol[:,1:2:end])')
plot(log.(sol[[1, 3],1:2:end])')
