


const neuts_max = 8
const k = 4
const n_neut_steps = neuts_max * k

const n_compartments = 3


const n_pop = 1.0
const n_inf_0 = 1.0 / 1000

const c_sus = 1
const c_exp = 2
const c_inf = 3


const beta = 1.0
const sigma = 1.0
const gamma = 0.5


const lambda = 0.05

const neut_levels = (collect(1:n_neut_steps)) / k
const protection_acquisition = 1 ./ (1 .+ exp.(-4 .* (neut_levels .- 2.0)))

const neut_jump_dist = Normal(5, 0.5)