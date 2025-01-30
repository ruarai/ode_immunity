

function get_sus(ode_solution, t, model_params)
    return Matrix(ode_solution(t)[ode_ix_sus(1:model_params.S), :]')
end

function get_inf(ode_solution, t, model_params)
    return ode_solution(t)[ode_ix_inf(model_params.S), :]
end

function get_inc(ode_solution, t, model_params; maintain_length = true)
    if maintain_length
        vcat([NaN], diff(ode_solution(t)[ode_ix_count(model_params.S), :]))
    else
        diff(ode_solution(t)[ode_ix_count(model_params.S), :])
    end
end

function get_results(ode_solution, t, model_params)
    return (
        sus = get_sus(ode_solution, t, model_params),
        inf = get_inf(ode_solution, t, model_params),
        inc = get_inc(ode_solution, t, model_params)
    )
end

function get_results_boosting(ode_solution, t, model_params)
    return (
        sus = Matrix(ode_solution(t)[ode_ix_boosting(c_sus, model_params.S, 1:model_params.S), :]'),
        inf = Matrix(ode_solution(t)[ode_ix_boosting(c_inf, model_params.S, 1:model_params.S), :]'),
        count = ode_solution(t)[model_params.S * 2 + 1, :]
    )
end

function get_summ_boosting(ode_solution, t, model_params)
    sus, inf, count = get_results_boosting(ode_solution, t, model_params)
    return (
        sus = sus, 
        inf = sum(inf, dims = 2), 
        inc = vcat([NaN], diff(count))
    )
end


function get_period(ode_sol, model_params, burn_in_days, n_days, Δt, ϵ)
    t = burn_in_days:Δt:n_days
    y = ode_sol(t)[1:(model_params.S + 1), :]'
    y_zeroed = copy(y)

    for j in axes(y, 1)
        y_zeroed[j, :] .= y[j, :] .- y[1, :]
    end
    sum_abs2 = vec(sum(abs2.(y_zeroed), dims = 2))

    t_repeats = t[findall(sum_abs2 .< ϵ)]
    n_repeats = length(t_repeats)
    t_repeat_groups = zeros(Int, length(t_repeats))
    t_repeat_groups[1] = 1

    for i in 2:n_repeats
        diff = t_repeats[i] - t_repeats[i - 1]
        if diff > 10.0
            t_repeat_groups[i] = t_repeat_groups[i - 1] + 1
        else
            t_repeat_groups[i] = t_repeat_groups[i - 1]
        end
    end

    t_repeats_grouped = zeros(maximum(t_repeat_groups))

    for i in eachindex(t_repeats_grouped)
        t_repeats_grouped[i] = mean(t_repeats[t_repeat_groups .== i])
    end

    period_samples = (t_repeats_grouped .- t_repeats_grouped[1]) ./ (0:(length(t_repeats_grouped) - 1))

    return (mean(period_samples[2:end]), std(period_samples[2:end]), length(period_samples[2:end]))
end


function get_periodic_attack_rate(ode_solution, model_params, burn_in_days, n_days, mean_period)
    if isnan(mean_period)
        return NaN
    end

    windows = [[i, i + mean_period] for i in burn_in_days:(n_days - mean_period)] 

    cumulative_attacks = [sum(ode_solution(t)[(model_params.S * 2 + 1):(model_params.S * 3),:], dims = 1) for t in windows]
    attack_mat = mapreduce(permutedims, hcat, cumulative_attacks)

    return mean(attack_mat[2, :] .- attack_mat[1, :])
end

function get_max_lyapunov(model_params, n_days, n_days_burn_in)
    u0 = zeros(Float64, n_compartments * model_params.S)
    u0[ode_ix(c_sus, 1, model_params.S)] = 1.0 - 0.01
    u0[ode_ix(c_inf, 1, model_params.S)] = 0.01
    
    
    dyn_system = CoupledODEs(ode_step_no_count!, u0, model_params)
    return lyapunov(dyn_system, n_days - n_days_burn_in, Ttr = n_days_burn_in, Δt = 100.0, d0 = 1e-9)
end

function get_seasonality_coordinates(ode_solution, t_seq, model_params)
    inc = get_inc(ode_solution, t_seq, model_params, maintain_length = false)

    t_mod = t_seq[2:end] .% 365
    theta = (t_mod ./ 365) * 2 * pi

    x = cos.(theta)
    y = sin.(theta)

    weighted_x = sum(x .* inc) / sum(inc)
    weighted_y = sum(y .* inc) / sum(inc)

    return (weighted_x, weighted_y)
end

