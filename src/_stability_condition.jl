

t = 0.01:0.01:100


X = Truncated(Normal(0.5, 0.5), 0, Inf)
plot(t, cdf.(X, t))

y_t = cdf.(X, t) ./ cdf.(X, 2 .* t)



plot(t, y_t)



means = 0.1:0.1:10
sds = 0.01:0.1:10


res = zeros(length(means), length(sds))
res2 = zeros(length(means), length(sds))

Threads.@threads for (i, μ) in collect(enumerate(means))
    for (j, σ) in enumerate(sds)

        X = Normal(μ, σ)
        y_t = cdf.(X, t) ./ cdf.(X, 2 .* t)

        res[i, j] = minimum(filter(!isnan,y_t))
        res2[i, j] = cdf.(X, 0)
    end
end


heatmap(sds, means, res)
heatmap(sds, means, res2)
heatmap(sds, means, res .> 1 / 4)
heatmap(sds, means, res2 .> 1 / 120)
Plots.abline!(2.20, -0.1, line=:dash)



X = Normal(3, 1)

t = 0.01:0.01:10

F_t = cdf.(X, t)
F_2t = cdf.(X, 2 .* t)

P_t = 1 .- cdf.(X, t)
P_2t = 1 .- cdf.(X, 2 .* t)

# P_(t)
plot(t, P_t)
plot!(t, P_2t)


plot(t, F_2t)
plot!(t, 4 .* F_t)


plot(t, F_t ./ F_2t)
plot!(t, F_t)


# 1 - P(2t) - 4 (1 - P(t)) ( > 0 )
plot(t, 1 .- P_2t .- 4 .* (1 .- P_t))


upp_t = (1 / sqrt(2 * π)) .* (1 ./ t) .* exp.((-t .^ 2) ./ 2)
low_t = (1 / sqrt(2 * π)) .* (1 ./ (t .^ 2 .+ 1)) .* exp.((-t .^ 2) ./ 2)

# P_(t), upper(t)
plot(t, P_t, ylim = (0, 1))
plot!(t, P_2t, ylim = (0, 1))
plot!(t, upp_t)
plot!(t, low_t)



means = 0.1:0.1:10

function dist(μ, σ)
    X = Normal(μ, σ)


    t = 0.01:0.01:10
    
    F_t = cdf.(X, t)
    F_2t = cdf.(X, 2 .* t)
    
    f_t = F_2t ./ F_t .- 4
    return abs(maximum(f_t))
end

μ = 3.0

using Optim
x0 = [0.0]

sd_opt = zeros(length(means))

for (i, μ) in enumerate(means)

    f(x) = dist(μ, exp(x[1]))
    res = optimize(f, -10, 100; rel_tol = 1e-16)
    sd_opt[i] = exp(Optim.minimizer(res)[1])
end

plot(means, sd_opt)

sd_opt ./ means





X = Normal(3, 0.5)

x = rand(X, 100)


x_t = zeros(length(x), length(t))
for i in axes(x_t, 1), j in axes(x_t, 2)
    x_t[i, j] = x[i] - 0.1 * t[j]
end

t_0 = zeros(length(x))

for i in axes(x_t, 1)
    t_0[i] = t[findfirst(x_t[i,:] .< 0)]
end

heatmap(x_t .< 0)

histogram(t_0)




γ = 0.5
μ = 0.1:1:35
σ = 1.0



lower = 2 .* (1 .+ μ .* γ) ./ (γ .^ 2 .* (σ .^ 2 .+ μ .^ 2)) .+ 1
upper = 2 .* μ .* γ .+ 3


plot(μ, lower)
plot!(μ, upper, ylim = (0, 2.5))

plot!([30], [1.5], seriestype = :scatter)
plot!([12], [1.1], seriestype = :scatter)