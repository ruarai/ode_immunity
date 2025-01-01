
struct model_parameters
    k::Int64
    S::Int64

    beta::Float64
    gamma::Float64

    eta::Float64

    C::Float64
    r::Float64
    rho::Float64

    wane_transition_rate::Float64

    b::Float64
    m::Float64

    c_jump_dist::Distribution

    c_levels::Vector{Float64}
    p_acq::Vector{Float64}

    M::Matrix{Float64}
    p_trans::Vector{Float64}
end



function make_model_parameters(;
    k,
    beta, gamma, 
    C, r,
    b, h,
    c_jump_dist,
    eta = 0.0,
    boosting = "independent"
)
    S = k + 1
    c_levels = collect(10 .^ (C .* (0:k) / k))
    
    p_acq = (c_levels .^ h) ./ (b ^ h .+ c_levels .^ h)

    p_trans = zeros(S)
    M = zeros(S, S)
    M[1, :] .= 1 # Null boosting matrix to make unit tests happy

    if boosting == "independent" # i.e., no boosting
        p_trans = build_immunity_transition_vector(S, c_levels, c_jump_dist)
    elseif boosting == "loglinear"
        M = build_immunity_matrix_boost_loglinear(S, c_levels, c_jump_dist)
    elseif boosting == "none" # Should be equivalent to independent, but allowing I stratification (and preventing decrease in strata)
        M = build_immunity_matrix_no_boost(S, c_levels, c_jump_dist)
    else
        throw(ArgumentError("Unknown boosting method specified"))
    end

    # Calculate rho from decay rate r
    rho = r / (C * log(10))

    wane_transition_rate = rho * k

    return model_parameters(
        k, S,
        beta, gamma, 
        eta,
        C, r, rho,
        wane_transition_rate,
        b, h,
        c_jump_dist,
        c_levels, p_acq,
        M, p_trans
    )    
end

# Create the post-infection immunity matrix M. Captures the probability of transitioning from
# strata j to strata i (bit backwards so matrix multiplication works later)
function build_immunity_matrix_boost_loglinear(N, c_levels, c_jump_dist)
    mat_immunity = zeros(N, N)

    log_c_levels = log10.(c_levels)

    for j in 1:N, i in j:N
        if j == N
            mat_immunity[i, j] = 1
        elseif i == j
            mat_immunity[i, j] = cdf(c_jump_dist, log_c_levels[i + 1] - log_c_levels[j])
        elseif i == N
            mat_immunity[i, j] = 1 - cdf(c_jump_dist, log_c_levels[i] - log_c_levels[j])
        else
            mat_immunity[i, j] = cdf(c_jump_dist, log_c_levels[i + 1] - log_c_levels[j]) - cdf(c_jump_dist, log_c_levels[i] - log_c_levels[j])
        end
    end

    return mat_immunity
end

function build_immunity_matrix_no_boost(N, c_levels, c_jump_dist)
    mat_immunity = zeros(N, N)

    log_c_levels = log10.(c_levels)

    for j in 1:N, i in j:N
        if j == N
            mat_immunity[i, j] = 1
        elseif i == j
            mat_immunity[i, j] = cdf(c_jump_dist, log_c_levels[i + 1])
        elseif i == N
            mat_immunity[i, j] = 1 - cdf(c_jump_dist, log_c_levels[i])
        else
            mat_immunity[i, j] = cdf(c_jump_dist, log_c_levels[i + 1]) - cdf(c_jump_dist, log_c_levels[i])
        end
    end

    return mat_immunity
end

# For the case without I stratification, i.e. no boosting
function build_immunity_transition_vector(N, c_levels, c_jump_dist)
    p_trans = zeros(N)
    log_c_levels = log10.(c_levels)

    for i in 1:N
        if i == N
            p_trans[i] = 1 - cdf(c_jump_dist, log_c_levels[i])
        elseif i == 1
            p_trans[i] = cdf(c_jump_dist, log_c_levels[i + 1])
        else
            p_trans[i] = cdf(c_jump_dist, log_c_levels[i + 1]) - cdf(c_jump_dist, log_c_levels[i])
        end
    end

    return p_trans
end