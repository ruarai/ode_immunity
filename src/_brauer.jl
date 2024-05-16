include("dependencies.jl")

using DataFrames


n_days = 365 * 4

c_jump_dist = Normal(1.0, 0.0)
lambda = 0.04
τ = 1.0 / lambda

γ = 0.5
R0 = 1.5
β = R0 * γ

t_days, S_t, I_t, R_t, dI_t = sim_dde(;
    R0 = R0,
    γ = γ,
    τ = c_jump_dist / lambda,
    dt = 0.01,
    n_days, 
    n_inf_0 = 0.001
)

plot(t_days, I_t)


function I_eq(γ, τ, β)
    return ((β .- γ) ./ (1 .+ τ .* γ)) ./ β
end
hline!([I_eq(γ, τ, β)])

plot(t_days, S_t)
hline!([γ/β])


# TODO: double check if there is an error in the minus one here?
eq_A(y, γ, τ, β) = (γ ./ y) .* sin.(τ .* y) .+ 1
eq_B(y, γ, τ, β) = (γ ./ y) .* β .* I_eq(γ, τ, β) .* (1 .- cos.(τ .* y)) .- y




using IntervalArithmetic, IntervalRootFinding


lambdas = 0.001:0.001:0.05
n_params = length(lambdas)
taus = 1 ./ lambdas

Beta = 0.2..1.5

betas = Array{Array{Float64}}(undef, n_params)
does_converge = zeros(Bool, n_params)

Threads.@threads for i in eachindex(taus)

    Y = (1e-4)..(2π / taus[i])

    γ = 0.5
    τ = taus[i]

    function root_fn((y, beta))
        return SVector(
            eq_A(y, γ, τ, beta) ^ 2 + eq_B(y, γ, τ, beta) ^ 2,
            0.0
        )
    end

    rts = roots(root_fn, Y × Beta, IntervalRootFinding.Krawczyk, 1e-8)

    if length(rts) > 1
        min_betas = zeros(length(rts))
        for j in eachindex(rts)
            min_betas[j] = sup(IntervalRootFinding.root_region(rts[j])[2])
        end
    
        betas[i] = min_betas
        does_converge[i] = true
    end
end

[length(unique(b)) for b in betas]

beta_firsts = [b[1] for b in betas]

plot!(beta_firsts[does_converge] / 0.5, lambdas[does_converge], xlim = (1.0, 3.0))

ys = 0:1e-4:0.5

plot(ys, eq_A(ys, γ, τ, min_beta), ylim = (-0.5, 0.5), xlim = (0, 0.5))
plot!(ys, eq_B(ys, γ, τ, min_beta))

vline!([min_y])










min_err = Inf
min_y = (0, 0)
min_beta = Inf

for i in eachindex(rts)
    root_box = IntervalRootFinding.root_region(rts[i])
    err = abs(sup(root_box[1]) - sup(root_box[2]))
    if err < min_err
        min_err = err
        min_y = (sup(root_box[1]), sup(root_box[2]))
        min_beta = sup(root_box[3])
    end
end

min_err
min_y
min_beta