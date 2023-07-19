
# Create the waning matrix. Captures (unscaled) waning when multiplied by vector of susceptibles.
function get_waning_matrix(n_neut_steps)
    mat_waning = zeros(n_neut_steps, n_neut_steps)

    for i in 1:n_neut_steps, j in 1:n_neut_steps
        if i == j - 1
            mat_waning[i, j] = 1
        elseif i == j && i != 1
            mat_waning[i, j] = -1
        end
    end

    return mat_waning
end

# Create the post-infection immunity matrix. Captures the probability of transitioning from
# strata i to strata j
function get_immunity_matrix(n_neut_steps, neut_levels, neut_jump_dist)
    mat_immunity = zeros(n_neut_steps, n_neut_steps)

    for i in 1:n_neut_steps, j in 1:n_neut_steps
        neuts_i = neut_levels[j]
        if i == 1
            neuts_j_m1 = 0.0
        else 
            neuts_j_m1 = neut_levels[i - 1]
        end
    
        neuts_j = neut_levels[i]
    
    
        if i == n_neut_steps
            # Edge case, captures upper bound of probability
            mat_immunity[i, j]  = 1 - cdf(neut_jump_dist, neuts_j_m1 - neuts_i)
        else
            mat_immunity[i, j] = cdf(neut_jump_dist, neuts_j - neuts_i) - cdf(neut_jump_dist, neuts_j_m1 - neuts_i)
        end
    end
    return mat_immunity
end

const B = get_waning_matrix(n_neut_steps)
const M = get_immunity_matrix(n_neut_steps, neut_levels, neut_jump_dist)
