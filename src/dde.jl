function sim_dde(; R0, γ, τ, dt, n_days, n_inf_0, return_each_step = true)
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

    if return_each_step
        return(0:dt:n_days, S_t, I_t, R_t, dI_t)
    else
        return (0:n_days,
                S_t[1:round(Int, 1 / dt):end],
                I_t[1:round(Int, 1 / dt):end],
                R_t[1:round(Int, 1 / dt):end],
                dI_t[1:round(Int, 1 / dt):end])
    end


end