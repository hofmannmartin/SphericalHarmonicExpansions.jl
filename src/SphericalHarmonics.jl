module SphericalHarmonics

using MultivariatePolynomials

include("sphericalHarmonics.jl")
export ylm

include("sphericalHarmonicsExpansion.jl")
export sphericalHarmonicsExpansion, SphericalHarmonicCoefficients

end # module
