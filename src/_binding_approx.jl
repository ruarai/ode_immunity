
include("dependencies.jl")

function ode_ab!(du, u, p, t)
    k1 = p[1]
    k2 = p[2]

    x, y, xy = u

    du[1] = -k1 * x * y + k2 * xy
    du[2] = -k1 * x * y + k2 * xy
    du[3] = k1 * x * y - k2 * xy
end

x_c1 = 10 .^ (-3:0.1:3)
x_c2 = 10 .^ (-3:0.5:3)
y_prop = zeros(length(x_c1), length(x_c2))

for i in eachindex(x_c1), j in eachindex(x_c2)
    
    u0 = zeros(Double64, 3)
    u0[1] = x_c1[i]
    u0[2] = x_c2[j]

    p = (1.0, 1.0)

    prob = ODEProblem{true, SciMLBase.FullSpecialize}(ODEFunction(ode_ab!), u0, (0.0, 50.0), p)

    sol = DifferentialEquations.solve(prob, Euler(), dt = 0.001, saveat = 0.01)

    y_prop[i, j] = sol[3,end] / (sol[2,end] + sol[3,end])
end


jldsave("data/paper/binding_approx.jld2"; x_c1, x_c2, y_prop)

plot(log10.(x_c1), y_prop)
vline!([log10.(1.0)])

heatmap(y_prop)
