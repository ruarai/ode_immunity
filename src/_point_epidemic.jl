


include("dependencies.jl")

k = 1
beta = 1.4
gamma = 1.0

lambda = 0.01

n_inf_0 = 0.001

# S0, I, R1 ... Rk
u0 = zeros(Double64, 2 + k)
u0[1] = 1 - n_inf_0
u0[2] = n_inf_0


function final_size(S0, R0)
    return S0 .+ (1 ./ R0) .* lambertw.(-R0 .* S0 .* exp.(-R0 .* S0))
end

final_sizes = [final_size(x, beta / gamma) for x in 0:0.001:1]


dt = 0.01

plot()


function step_fixed!(du, u, params, dt)

    S = u[1]
    I = u[2]
    Rn = u[3:(2+k)]

    β, γ, λ = params

    epi_size = 0
    if S > 1.2 / (β / γ)
        epi_size = final_sizes[round(Int, S * 1000) + 1]
    end

    dS = - epi_size + Rn[1] * λ * dt
    dI = 0

    dRk = epi_size - Rn[k] * λ * dt
    dRn = Rn[2:k] * λ * dt - Rn[1:(k - 1)] * λ * dt


    du[1] = dS
    du[2] = dI
    du[3:(2 + k - 1)] = dRn
    du[2 + k] = dRk
end

du = zeros(Double64, length(u0))

params = (beta, gamma, lambda)


n_days = 500
u = copy(u0)
u_t = zeros(n_days, length(u0))

for t in 1:n_days
    u_t[t,:] = u
    for s in 1:100
        step_fixed!(du, u, params, dt)
        u = u + du
    end
end

plot(u_t)
