module SphericalHarmonics

using MultivariatePolynomials

include("sphericalHarmonics.jl")
export ylm

include("sphericalHarmonicsExpansion.jl")
export sphericalHarmonicsExpansion, SphericalHarmonicCoefficients

include("fasteval.jl")
export @fasteval

include("sphericalQuadrature.jl")
export sphericalQuadrature

end # module
