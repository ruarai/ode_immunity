
# Create the waning matrix B. Captures (unscaled) waning when multiplied by vector of susceptibles.
function build_waning_matrix(params)
    N = params.N

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
function build_immunity_matrix(params)
    N = params.N
    c_levels = params.c_levels
    c_jump_dist = params.c_jump_dist

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
