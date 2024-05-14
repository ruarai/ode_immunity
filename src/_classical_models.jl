
include("dependencies.jl")

using JLD2

k = 16
R = 1.5
gamma = 0.25
beta = R * gamma
lambda = 0.003
b = 0.25
m = 40
c_jump_dist = Normal(0.8, 0.05)

model_params = make_model_parameters(
    k = k, beta = beta, gamma = gamma, lambda = lambda,
    b = b, m = m, c_jump_dist = c_jump_dist; boosting = false
)

ode_sparsity = ode_get_sparsity(model_params)

n_inf_0 = 0.0001
n_days = 365*20

ode_solution = @time ode_solve(model_params, n_days, n_inf_0, ode_sparsity)

c_levels = model_params.c_levels


plot(sum(sol_I, dims = 2)[1:1000])


model_params_classical = (
    beta = model_params.beta,
    gamma = model_params.gamma,
    eta = 0.01
)

function ode_step_si!(du, u, model_params, t)
    β = model_params.beta

    flow_inf = β * u[1] * u[2]
    du[1] = -flow_inf
    du[2] = +flow_inf
end

u0 = [1 - n_inf_0, n_inf_0]
tspan = (0.0, 1.0 * n_days)
prob = ODEProblem{true, SciMLBase.FullSpecialize}(ode_step_si!, u0, tspan, model_params_classical)

sol_si = DifferentialEquations.solve(prob, Euler(), dt = 0.01, saveat = 1);

plot(sol_si)



function ode_step_sir!(du, u, model_params, t)
    β = model_params.beta
    γ = model_params.gamma

    flow_inf = β * u[1] * u[2]
    flow_rec = γ * u[2]
    du[1] = -flow_inf
    du[2] = +flow_inf - flow_rec
    du[3] = +flow_rec
end


u0 = [1 - n_inf_0, n_inf_0, 0]
tspan = (0.0, 1.0 * n_days)
prob = ODEProblem{true, SciMLBase.FullSpecialize}(ode_step_sir!, u0, tspan, model_params_classical)

sol_sir = DifferentialEquations.solve(prob, Euler(), dt = 0.01, saveat = 1);



function ode_step_sirs!(du, u, model_params, t)
    β = model_params.beta
    γ = model_params.gamma
    η = model_params.eta

    flow_inf = β * u[1] * u[2]
    flow_rec = γ * u[2]
    flow_waning = η * u[3]
    du[1] = flow_waning - flow_inf
    du[2] = flow_inf - flow_rec
    du[3] = flow_rec - flow_waning
end


u0 = [1 - n_inf_0, n_inf_0, 0]
tspan = (0.0, 1.0 * n_days)
prob = ODEProblem{true, SciMLBase.FullSpecialize}(ode_step_sirs!, u0, tspan, model_params_classical)

sol_sirs = DifferentialEquations.solve(prob, Euler(), dt = 0.01, saveat = 1);



function sim_dde(; R0, γ, τ, dt, n_days, n_inf_0)
    t_steps = 0:dt:n_days
    n_t_steps = length(t_steps)

    β = R0 * γ

    S_t = zeros(n_t_steps)
    I_t = zeros(n_t_steps)
    R_t = zeros(n_t_steps)
    dI_t = zeros(n_t_steps)

    S_t[1] = 1 - n_inf_0
    I_t[1] = n_inf_0

    P_t = 1 .- cdf.(τ, t_steps)
    ix_P_t_neglible = findfirst(P_t .< 1e-8)
    if isnothing(ix_P_t_neglible)
        ix_P_t_neglible = length(P_t)
    end

    println(ix_P_t_neglible)
    
    for t in 1:(n_t_steps - 1)
        I = I_t[t]
    
        recovered = 0
        if t > 1
            for x in 1:min(t - 1, ix_P_t_neglible)
                recovered += @views I_t[t - x] * P_t[x + 1] * dt
            end
        end

        R = γ * recovered

    
        S = 1 - I - R

    
        dI = -γ * I + β * I * S
    
        S_t[t + 1] = S - dI * dt
        I_t[t + 1] = I + dI * dt
        R_t[t + 1] = 1 -  S_t[t + 1] - I_t[t + 1]
        dI_t[t + 1] = dI 

    end

    return (0:n_days,
            S_t[1:round(Int, 1 / dt):end],
            I_t[1:round(Int, 1 / dt):end],
            R_t[1:round(Int, 1 / dt):end],
            dI_t[1:round(Int, 1 / dt):end])
end

c_jump_dist = Normal(0.8 - 0.25, 0.05)
lambda = 0.003

τ = c_jump_dist / lambda

t_days, S_t, I_t, R_t, dI_t = sim_dde(;
    R0 = model_params.beta / model_params.gamma ,
    γ = model_params.gamma,
    τ = τ,
    dt = 0.01,
    n_days, 
    n_inf_0 = n_inf_0
)

plot(I_t)
plot!(sol_I)

plot(sol_si[2,:])
plot(sol_sir[2,:])
plot!(sol_sirs[2,:])


I_structured = zeros(n_days)
I_dde = zeros(n_days)
I_SI = zeros(n_days)
I_SIR = zeros(n_days)
I_SIRS = zeros(n_days)

for d in 1:n_days
    I_structured[d] = sum(ode_solution(d)[ode_ix(c_inf, 1:model_params.S, model_params.S)])
    I_dde[d] = I_t[d]
    I_SI[d] = sol_si(d)[2]
    I_SIR[d] = sol_sir(d)[2]
    I_SIRS[d] = sol_sirs(d)[2]
end


plot(I_structured)
plot!(I_dde)
#plot!(I_SI)
plot!(I_SIR)
plot!(I_SIRS)



jldsave("data/paper/classical.jld2"; I_structured, I_dde, I_SI, I_SIR, I_SIRS)


