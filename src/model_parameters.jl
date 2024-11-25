
struct model_parameters
    k::Int64
    S::Int64

    beta::Float64
    gamma::Float64

    eta::Float64

    C::Float64
    rho::Float64

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

    beta, gamma, 

    
    C, rho,

    b, h,

    c_jump_dist,

    eta = 0.0,
    boosting = "none"
)
    S = k + 1
    c_levels = collect(10 .^ (C .* (0:k) / k))
    
    p_acq = (c_levels .^ h) ./ (b ^ h .+ c_levels .^ h)

    B = build_waning_matrix(S)
    if boosting == "linear"
        M = build_immunity_matrix_boost(S, c_levels, c_jump_dist)
    elseif boosting == "loglinear"
        M = build_immunity_matrix_boost_loglinear(S, c_levels, c_jump_dist)
    elseif boosting == "none"
        M = build_immunity_matrix_no_boost(S, c_levels, c_jump_dist)
    else
        throw(ArgumentError("Unknown boosting method specified"))
    end

    wane_transition_rate = rho * k

    return model_parameters(
        k, S,

        beta, gamma, 

        eta,
        
        C, rho,

        wane_transition_rate,

        b, h,
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
# strata j to strata i (bit backwards so matrix multiplication works later)
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

function build_immunity_matrix_boost_loglinear(N, c_levels, c_jump_dist)
    mat_immunity = zeros(N, N)

    for j in 1:N, i in j:N
        if j == N
            mat_immunity[i, j] = 1
        elseif i == j
            mat_immunity[i, j] = cdf(c_jump_dist, log2(c_levels[i + 1]) - log2(c_levels[j]))
        elseif i == N
            mat_immunity[i, j] = 1 - cdf(c_jump_dist, log2(c_levels[i]) - log2(c_levels[j]))
        else
            mat_immunity[i, j] = cdf(c_jump_dist, log2(c_levels[i + 1]) - log2(c_levels[j])) - cdf(c_jump_dist, log2(c_levels[i]) - log2(c_levels[j]))
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