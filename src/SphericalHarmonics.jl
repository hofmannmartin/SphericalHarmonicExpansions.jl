__precompile__()
module SphericalHarmonics

using Reexport
using MultivariatePolynomials
@reexport using TypedPolynomials

include("sphericalHarmonics.jl")
export ylm

include("sphericalHarmonicsExpansion.jl")
export sphericalHarmonicsExpansion, SphericalHarmonicCoefficients

include("fastfunc.jl")
export @fastfunc

include("sphericalQuadrature.jl")
export sphericalQuadrature

include("translation.jl")
export addRl, rlm

end # module
