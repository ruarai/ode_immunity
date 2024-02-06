include("dependencies.jl")



function ode_step_sirws!(du, u, sirws_params, t)
   
    S = u[1]
    I = u[2]
    R = u[3]
    W = u[4]

    β, γ, K, v = sirws_params 

    du[1] = - β * S * I + 2 * K * W
    du[2] = β * S * I - γ * I
    du[3] = γ * I - 2 * K * R + v * β * I * W
    du[4] = 2 * K * R - v * β * I * W - 2 * K * W
end


n_days = 10000
u0 = zeros(4)
u0[1] = 0.99
u0[2] = 0.01

sirws_params = (1.1, 0.5, 0.01, 10.0)

ode_step_fn = ODEFunction(ode_step_sirws!)

tspan = (0.0, 1.0 * n_days)
prob = ODEProblem{true, SciMLBase.FullSpecialize}(ode_step_fn, u0, tspan, sirws_params)

sol = DifferentialEquations.solve(prob, Euler(), dt = 0.01, saveat = 1);



plot(sol)


# Without "boosting"
function ode_step_sirws!(du, u, sirws_params, t)
   
    S = u[1]
    I = u[2]
    R = u[3]
    W = u[4]

    β, γ, a, b = sirws_params 

    du[1] = - β * S * I + b * W
    du[2] = β * S * I - γ * I
    du[3] = γ * I - a * R
    du[4] = a * R - b * W
end


n_days = 40000
u0 = zeros(4)
u0[1] = 0.99
u0[2] = 0.01

sirws_params = (1.1, 0.5, 0.001, 0.1)

ode_step_fn = ODEFunction(ode_step_sirws!)

tspan = (0.0, 1.0 * n_days)
prob = ODEProblem{true, SciMLBase.FullSpecialize}(ode_step_fn, u0, tspan, sirws_params)

sol = DifferentialEquations.solve(prob, Euler(), dt = 0.001, saveat = 1);



plot(sol)