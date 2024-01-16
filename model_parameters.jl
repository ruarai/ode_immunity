
struct model_parameters
    c_max::Int32
    k::Int32

    beta::Float64
    sigma::Float64
    gamma::Float64

    lambda::Float64

    b::Float64
    m::Float64

    c_jump_dist::Distribution

    N::Int64

    c_levels::Vector{Float64}
    p_acq::Vector{Float64}

    B::Matrix{Float64}
    M::Matrix{Float64}
end



function make_model_parameters(;
    c_max, k,

    beta, sigma, gamma, lambda,

    b, m,

    c_jump_dist
)
    N = c_max * k
    c_levels = collect(1:N) ./ k
    p_acq = 1 ./ (1 .+ exp.(-m .* (c_levels .- b)))

    B = build_waning_matrix(N)
    M = build_immunity_matrix(N, c_levels, c_jump_dist)

    return model_parameters(
        c_max, k,

        beta, sigma, gamma, lambda,

        b, m,
        c_jump_dist,

        N,
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
function build_immunity_matrix(N, c_levels, c_jump_dist)
    mat_immunity = zeros(N, N)

    for i in 1:N, j in 1:N
    
        # Edge cases for upper and lower bounds
        if i == N
            mat_immunity[i, j]  = 1 - cdf(c_jump_dist, c_levels[i - 1] - c_levels[j])
        elseif i == 1
            mat_immunity[i, j]  = cdf(c_jump_dist, c_levels[i] - c_levels[j])
        else
            mat_immunity[i, j] = cdf(c_jump_dist, c_levels[i] - c_levels[j]) - cdf(c_jump_dist, c_levels[i - 1] - c_levels[j])
        end
    end
    return mat_immunity
end