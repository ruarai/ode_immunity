
struct model_parameters
    k::Int64

    beta::Float64
    gamma::Float64

    lambda::Float64

    wane_transition_rate::Float64

    b::Float64
    m::Float64

    c_jump_dist::Distribution

    c_levels::Vector{Float64}
    p_acq::Vector{Float64}

    B::Matrix{Float64}
    M::Matrix{Float64}
end



function make_model_parameters(;
    k,

    beta, gamma, lambda,

    b, m,

    c_jump_dist,

    boosting = true
)
    c_levels = collect((0 : (k - 1)) / k)
    
    p_acq = 1 ./ (1 .+ exp.(-m .* (c_levels .- b)))

    B = build_waning_matrix(k)
    if boosting
        M = build_immunity_matrix_boost(k, c_levels, c_jump_dist)
    else
        M = build_immunity_matrix_no_boost(k, c_levels, c_jump_dist)
    end

    wane_transition_rate = lambda * (k - 1)

    return model_parameters(
        k,

        beta, gamma, lambda,

        wane_transition_rate,

        b, m,
        c_jump_dist,

        c_levels, p_acq,

        B, M
    )    
end


# Create the waning matrix B. Captures (unscaled) waning when multiplied by vector of susceptibles.
function build_waning_matrix(N)
    mat_waning = zeros(N, N)

    for i in 1:N, j in 1:N
        if i == j - 1
            mat_waning[i, j] = 1
        elseif i == j && i != 1
            mat_waning[i, j] = -1
        end
    end

    return mat_waning
end

# Create the post-infection immunity matrix M. Captures the probability of transitioning from
# strata j to strata i (bit backwards so matrix multiplication works)
function build_immunity_matrix_boost(N, c_levels, c_jump_dist)
    mat_immunity = zeros(N, N)

    for j in 1:N, i in j:N
        if j == N
            mat_immunity[i, j] = 1
        elseif i == j
            mat_immunity[i, j] = cdf(c_jump_dist, c_levels[i + 1] - c_levels[j])
        elseif i == N
            mat_immunity[i, j] = 1 - cdf(c_jump_dist, c_levels[i] - c_levels[j])
        else
            mat_immunity[i, j] = cdf(c_jump_dist, c_levels[i + 1] - c_levels[j]) - cdf(c_jump_dist, c_levels[i] - c_levels[j])
        end
    end

    return mat_immunity
end


function build_immunity_matrix_no_boost(N, c_levels, c_jump_dist)
    mat_immunity = zeros(N, N)

    for j in 1:N, i in j:N
        if j == N
            mat_immunity[i, j] = 1
        elseif i == j
            mat_immunity[i, j] = cdf(c_jump_dist, c_levels[i + 1])
        elseif i == N
            mat_immunity[i, j] = 1 - cdf(c_jump_dist, c_levels[i])
        else
            mat_immunity[i, j] = cdf(c_jump_dist, c_levels[i + 1]) - cdf(c_jump_dist, c_levels[i])
        end
    end

    return mat_immunity
end