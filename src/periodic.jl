function get_period(y)
    period_lens = 10:(length(y) ÷ 2)
    err_period = zeros(length(period_lens))

    log_y = log.(y)

    min_err = Inf
    min_err_index = -1
    for i in eachindex(period_lens)
        err_p = 0.0

        diff = abs.(log_y[1:(end - period_lens[i]), :]' .- log_y[(period_lens[i] + 1):(end), :]')

        for j in eachindex(diff)
            if !isnan(diff[j]) && !isinf(diff[j])
                err_p += diff[j]
            end
        end

        err_period[i] = err_p

        if err_p < min_err - 20
            min_err = err_p
            min_err_index = i
        end
    end

    
    return period_lens[min_err_index]
end

function get_periodic_attack_rate(incidence, period_len)
    windows = [(i, i + period_len - 1) for i in 1:(length(incidence) - period_len)]
    sums = [sum(incidence[w[1]:w[2]]) for w in windows]

    if length(sums) >= 1
        return median(sums)
    end
    return NaN
end