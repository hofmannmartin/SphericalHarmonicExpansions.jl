module SphericalHarmonics

using LinearAlgebra, Reexport, HDF5, GeneralizedGenerated
@reexport using MultivariatePolynomials
@reexport using TypedPolynomials

include("sphericalHarmonic.jl")
export ylm, rlylm

include("sphericalHarmonicsExpansion.jl")
export sphericalHarmonicsExpansion, SphericalHarmonicCoefficients, solidHarmonicsExpansion,
        solid!, spherical!

include("fastfunc.jl")
export @fastfunc, fastfunc

include("sphericalQuadrature.jl")

include("translation.jl")
export translation

end # module
