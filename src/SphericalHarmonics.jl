module SphericalHarmonics

using MultivariatePolynomials

include("sphericalHarmonics.jl")
export ylm

include("sphericalHarmonicsExpansion.jl")
export sphericalHarmonicsExpansion, SphericalHarmonicCoefficients

include("fastfunc.jl")
export @fastfunc

include("sphericalQuadrature.jl")
export sphericalQuadrature

end # module
