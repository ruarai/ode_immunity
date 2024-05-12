using DifferentialEquations


include("dependencies.jl")


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


n_days = 365*20

c_jump_dist = Normal(1.0, 0.05)
lambda = 0.03

τ = c_jump_dist / lambda

t_days, S_t, I_t, R_t, dI_t = sim_dde(;
    R0 = 1.5,
    γ = 0.5,
    τ = τ,
    dt = 0.01,
    n_days, 
    n_inf_0 = 0.001
)

plot(t_days, I_t)
plot(t_days, S_t)
plot(t_days, R_t)
plot(t_days, dI_t)

plot(t_days, S_t .+ I_t .+ R_t)


lambdas = 0.01:0.001:0.05

results = zeros(length(lambdas), n_days + 1)

Threads.@threads for (i, lambda) in collect(enumerate(lambdas))
    println(lambda)

    τ = c_jump_dist / lambda
    t_days, S_t, I_t, R_t, dI_t = sim_dde(;
        R0 = 1.5, 
        γ = 0.5,
        τ = τ,
        dt = 0.01,
        n_days, 
        n_inf_0 = 0.01
    )
    results[i,:] = I_t
end

heatmap(results)
plot(results')
plot(results[findfirst(collect(lambdas) .== 0.0375),:])



periods = zeros(length(lambdas))
changes = zeros(length(lambdas))
sum_change = zeros(length(lambdas))

for i in eachindex(lambdas)
    res_subset = results[i, 2000:end]
    res_diff = diff(res_subset)
    # Adjust for half a day offset here?
    res_match = (res_diff[1:(end - 1)] .> 0) .& (res_diff[2:end] .<= 0)
    
    
    res_ix = findall(res_match)
    
    res_ix_offset = res_ix[2:end] .- res_ix[1]
    
    periods[i] = mean(res_ix_offset ./ (1:length(res_ix_offset)))

    a = res_subset[res_ix]

    changes[i] = mean(abs2.(a .- mean(a)))
    sum_change[i] = sum(res_diff)
end

plot(lambdas, periods)
plot(lambdas, changes)
plot(lambdas, abs.(sum_change) .< 1e-6)



i = 5
res_diff = diff(results[i, 2000:end])
res_match = (res_diff[1:(end - 1)] .> 0) .& (res_diff[2:end] .<= 0)

plot(results[i,2000:end])
plot!(res_match ./ 20)

res_ix = findall(res_match)

a = results[i, 2000:end][res_ix]
plot(a)

mean(abs.(a .- mean(a)))

res_ix_offset = res_ix[2:end] .- res_ix[1]

periods[i] = mean(res_ix_offset ./ (1:length(res_ix_offset)))