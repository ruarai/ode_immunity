
include("dependencies.jl")

using JLD2

k = 32
C = 8.0
R = 1.5
gamma = 0.25
beta = R * gamma
rho = 0.003
b = 2^3
h = 8
c_jump_dist = Normal(2^6, 2^5)

model_params = make_model_parameters(
    k = k, beta = beta, gamma = gamma, C = C, rho = rho,
    b = b, h = h, c_jump_dist = c_jump_dist; boosting = "none", eta = 0.15
)

ode_sparsity = ode_get_sparsity(model_params)

n_inf_0 = 0.0001
n_days = 365*100

ode_step_fn = ODEFunction(ode_step!; jac_prototype = float.(ode_sparsity))

u0 = zeros(Float64, n_compartments * model_params.S)
u0[ode_ix(c_sus, 1, model_params.S)] = 1.0 - n_inf_0
u0[ode_ix(c_inf, 1, model_params.S)] = n_inf_0

tspan = (0.0, 1.0 * n_days)
prob = ODEProblem{true, SciMLBase.FullSpecialize}(ode_step_fn, u0, tspan, model_params)




x_freq = 0.1:0.1:1.0
x_offset = 0.0:0.1:0.9



x_vals = vec([(freq = x1, offset = x2) for x1 in x_freq, x2 in x_offset])

y_inf = zeros(length(x_vals), n_days)
y_count = zeros(length(x_vals), n_days)

@showprogress Threads.@threads for i in eachindex(x_vals)
    vacc_gap = 365 * x_vals[i].freq
    t_vacc = (0:vacc_gap:n_days) .+ (x_vals[i].offset * vacc_gap)

    yearly_vacc_proportion = 0.4
    vacc_per_gap = (vacc_gap / 365) * yearly_vacc_proportion

    function cb_condition(u, t, integrator)
        return in(t, t_vacc)
    end

    function cb_affect!(integrator)
        to_vacc = 0.0
        vacc_prop = vacc_per_gap

        for i in 1:model_params.S
            to_vacc += integrator.u[ode_ix(c_sus, i, model_params.S)] * vacc_prop
            integrator.u[ode_ix(c_sus, i, model_params.S)] -= integrator.u[ode_ix(c_sus, i, model_params.S)] * vacc_prop

            to_vacc += integrator.u[ode_ix(c_inf, i, model_params.S)] * vacc_prop
            integrator.u[ode_ix(c_inf, i, model_params.S)] -= integrator.u[ode_ix(c_inf, i, model_params.S)] * vacc_prop
        end
        

        integrator.u[ode_ix(c_sus, model_params.S, model_params.S)] += to_vacc
    end

    cb = DiscreteCallback(cb_condition, cb_affect!)


    ode_solution = DifferentialEquations.solve(
        prob, Euler(), 
        callback = cb, tstops = t_vacc,
        dt = 0.005, saveat = 1
    );


    y_inf[i, :] = vec(sum(ode_solution(1:n_days)[ode_ix(c_inf, 1:model_params.S, model_params.S), :], dims = 1))
    y_count[i, :] = vec(sum(ode_solution(1:n_days)[ode_ix(c_count, 1:model_params.S, model_params.S), :], dims = 1))
end


plot(y_inf[:,:]')
plot(y_inf[2:end,:]')


y_max = maximum(y_inf[:,(365*80):end], dims = 2)

heatmap(reshape(y_max, length(x_freq), length(x_offset)))

plot(x_freq, reshape(y_max, length(x_freq), length(x_offset)))


plot(x_freq, maximum(y_inf[:,(365*80):end], dims = 2) / maximum(y_inf[end,(365*80):end]))

heatmap(y_diff)