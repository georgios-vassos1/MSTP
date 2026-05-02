using Distributions, LinearAlgebra, BlockDiagonals, Random

# Build a correlation block of size n with a single off-diagonal value
function corrblock(n::Int64, lower::Float64=0.6, upper::Float64=0.8, by::Float64=0.1)
    value = rand(lower:by:upper)
    [i == j ? 1.0 : value for i in 1:n, j in 1:n]
end

# Build a block-diagonal covariance matrix with optional cross-block correlation
function gen_cov_mat(n_blocks::Int64, block_size::Int64, cross_corr::Float64=0.0)
    blocks    = [corrblock(block_size) for _ in 1:n_blocks]
    covmat    = Matrix(BlockDiagonal(blocks))
    covmat[covmat .== 0.0] .= cross_corr
    covmat
end

# Clamp CDF values to (ε, 1−ε) so that Poisson quantile never returns Inf.
# Gaussian tails beyond ≈ ±8σ map to floating-point 0/1, causing infinite
# Poisson quantiles and LP infeasibility in the SDDP subproblems.
const _PCLAMP = 1e-10

# Sample N correlated Poisson scenario vectors — returns Vector{Vector{Int64}}
# Used for SDDP training: each element is one scenario (nOrigins + nDestinations values)
function sample_scenarios(N::Int64, lambda::Vector{Float64}, corrmat::Matrix{Float64})
    [vec(quantile.(Poisson.(lambda), clamp.(cdf.(Normal(), rand(MvNormal(corrmat), 1)), _PCLAMP, 1.0 - _PCLAMP))) for _ in 1:N]
end

# Sample N correlated Poisson scenarios — returns Matrix (N × dim)
# Used for bulk/OOB sampling
function sample_scenarios_mat(N::Int64, lambda::Vector{Float64}, corrmat::Matrix{Float64})
    quantile.(Poisson.(lambda), clamp.(cdf.(Normal(), rand(MvNormal(corrmat), N)), _PCLAMP, 1.0 - _PCLAMP))'
end
