

function get_sus(ode_solution, t, model_params)
    return Matrix(ode_solution(t)[ode_ix_sus(1:model_params.S), :]')
end

function get_inf(ode_solution, t, model_params)
    return ode_solution(t)[ode_ix_inf(model_params.S), :]
end

function get_inc(ode_solution, t, model_params)
    vcat([0], diff(ode_solution(t)[ode_ix_count(model_params.S), :]))
end

function get_results(ode_solution, t, model_params)
    return (
        sus = get_sus(ode_solution, t, model_params),
        inf = get_inf(ode_solution, t, model_params),
        inc = get_inc(ode_solution, t, model_params)
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