module SphericalHarmonics

using MultivariatePolynomials

include("sphericalHarmonics.jl")
include("sphericalHarmonicsExpansion.jl")
export ylm, sphericalHarmonicsExpansion

include("fasteval.jl")
export @fasteval

end # module
