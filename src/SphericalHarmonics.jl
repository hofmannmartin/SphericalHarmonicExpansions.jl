module SphericalHarmonics

using MultivariatePolynomials

include("sphericalHarmonics.jl")
include("sphericalHarmonicsExpansion.jl")
export ylm, sphericalHarmonicsExpansion

end # module
