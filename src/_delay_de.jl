using DifferentialEquations



n_inf_0 = 0.01

u0 = zeros(3)



function sim_dde(; R0, γ, τ, dt, n_days, n_inf_0)
    t_steps = 0:dt:n_days
    n_t_steps = length(t_steps)

    β = R0 * γ
    I_t = zeros(n_t_steps)
    I_t[1] = n_inf_0
    
    Pt = (1 .- cdf.(τ, reverse(t_steps)))
    
    for t in 1:(n_t_steps - 1)
        I = I_t[t]
    
        recovered = 0
        if t > 1
            P = Pt[(end - t + 1:(end - 1))]
            recovered = sum(P .* I_t[1:(t - 1)] .* dt)
        end
    
        S = 1 - I - γ * recovered
    
        dI = -γ * I + β * I * S
    
        I_t[t + 1] = I + dI * dt
    end

    return (1:n_days, I_t[1:round(Int, 1 / dt):(end - 1)])
end


n_days = 500

t_days, I_t = sim_dde(;
    R0 = 1.5, 
    γ = 0.5,
    τ = Normal(25, 1.0),
    dt = 0.01,
    n_days, 
    n_inf_0 = 0.01
)


plot(t_days, I_t)


means = 5:5:30

results = zeros(length(means), length(I_t))

Threads.@threads for (i, μ) in collect(enumerate(means))
    println(μ)
    _, I_t = sim_dde(;
        R0 = 1.5, 
        γ = 0.5,
        τ = Normal(μ, 1.0),
        dt = 0.01,
        n_days, 
        n_inf_0 = 0.01
    )
    results[i,:] = I_t
end

heatmap(t_days, means, log.(results))

plot(results')


