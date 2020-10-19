module SphericalHarmonics

using LinearAlgebra, Reexport, HDF5, GeneralizedGenerated
@reexport using MultivariatePolynomials
@reexport using TypedPolynomials

include("sphericalHarmonic.jl")
export ylm, rlm

include("sphericalHarmonicsExpansion.jl")
export sphericalHarmonicsExpansion, SphericalHarmonicCoefficients, solidHarmonicsExpansion,
        solid!, spherical!

include("fastfunc.jl")
export @fastfunc, fastfunc

include("sphericalQuadrature.jl")
export sphericalQuadrature, errorSphericalQuadrature

include("translation.jl")
export translation, errorTranslation

end # module
