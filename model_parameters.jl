
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

    N::Int32

    c_levels::Vector{Float64}
    p_acq::Vector{Float64}
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

    return model_parameters(
        c_max, k,

        beta, sigma, gamma, lambda,

        b, m,
        c_jump_dist,

        N,
        c_levels, p_acq
    )    
end