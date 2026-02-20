include("dependencies.jl")


function sir_ode!(du, u, p, t)
    S, I, R = u

    beta_t = p.beta * (1 + p.eta * cos(2π * t / 365.0))

    du[1] = -beta_t * S * I + p.sigma * R
    du[2] = beta_t * S * I - p.gamma * I
    du[3] = p.gamma * I - p.sigma * R
end


n_days_burn_in = 100 * 365
n_days = n_days_burn_in + 100 * 365


x_eta = 0.00:0.05:0.5
length(x_eta)

sigma_step = 0.001
x_sigma = sigma_step:sigma_step:0.01
length(x_sigma)

x_vals = vec([(eta = x1, sigma = x2) for x1 in x_eta, x2 in x_sigma])

y_period = zeros(3, length(x_vals))
y_inf_summary = zeros(4, length(x_vals))


u0 = [0.99, 0.01, 0.0]

@showprogress Threads.@threads for i in eachindex(x_vals)
    p = (beta = baseline_beta, gamma = baseline_gamma, sigma = x_vals[i].sigma, eta = x_vals[i].eta)
    tspan = (0.0, n_days) 

    prob = ODEProblem(sir_ode!, u0, tspan, p)

    sol = solve(
        prob,
        Rodas5P(),

        dtmax = 1.0,
        reltol = 1e-10, abstol = 1e-10,

        saveat = n_days_burn_in:0.25:n_days
    )

    y_period[:, i] .= get_period(sol, n_days_burn_in, n_days, 0.25, periodic_ϵ, 1:3)

    inf = [sol(t)[2] for t in n_days_burn_in:n_days]

    y_inf_summary[1, i] = minimum(inf)
    y_inf_summary[2, i] = maximum(inf)
    y_inf_summary[3, i] = mean(inf)
    y_inf_summary[4, i] = testchaos01(NaNMath.log10.(inf[1:80:end]))
end


y_per = y_period[1, :]

is_periodic = (y_period[1, :] .% 365 .< 1) .| (y_period[1, :] .% 365 .> 364)

y_per[.!is_periodic] .= NaN

y = reshape(y_per, length(x_eta), length(x_sigma))

heatmap(x_eta, x_sigma, min.(4, y' ./ 365))


y = reshape(y_inf_summary[1, :], length(x_eta), length(x_sigma))

heatmap(x_eta, x_sigma, max.(-10,log10.(y')))

