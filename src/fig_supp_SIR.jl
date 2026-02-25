include("dependencies.jl")


function sir_ode!(du, u, p, t)
    S, I, R = u

    beta_t = p.beta * (1 + p.eta * cos(2π * t / 365.0))

    du[1] = -beta_t * S * I + p.sigma * R
    du[2] = beta_t * S * I - p.gamma * I
    du[3] = p.gamma * I - p.sigma * R
end


n_days_burn_in = 1000 * 365
n_days = n_days_burn_in + 250 * 365

u0 = [0.99, 0.01, 0.0]
p = (beta = baseline_beta, gamma = baseline_gamma, sigma = 0.01, eta = 0.25)
tspan = (0.0, n_days) 



prob = ODEProblem(sir_ode!, u0, tspan, p)

sol = solve(
    prob,
    Rodas5P(),

    dtmax = 1.0,
    reltol = 1e-10, abstol = 1e-10,

    saveat = n_days_burn_in:0.25:n_days
)


get_period(
    sol, n_days_burn_in, n_days, 
    0.25, 1e-9,
    1:3
)

plot([sol(t)[3] for t in (n_days_burn_in + 1):(n_days_burn_in + 365 * 25)])

scatter(repeat(1:365, outer = 250), [sol(t)[2] for t in (n_days_burn_in + 1):(n_days_burn_in + 365 * 250)])