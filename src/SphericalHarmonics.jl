__precompile__()
module SphericalHarmonics

using Reexport, HDF5

@reexport using MultivariatePolynomials

if !isdir(Pkg.dir("TypedPolynomials"))
    println("Installing TypedPolynomials...")
    Pkg.clone("https://github.com/rdeits/TypedPolynomials.jl.git")
end
@reexport using TypedPolynomials

include("sphericalHarmonic.jl")
export ylm, rlm

include("sphericalHarmonicsExpansion.jl")
export sphericalHarmonicsExpansion, SphericalHarmonicCoefficients, solidHarmonicsExpansion,
        solid, spherical

include("fastfunc.jl")
export @fastfunc

include("sphericalQuadrature.jl")
export sphericalQuadrature, errorSphericalQuadrature

include("translation.jl")
export translation, translateRlm, errorTranslation

end # module
