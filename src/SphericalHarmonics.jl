__precompile__()
module SphericalHarmonics

using Reexport

using MultivariatePolynomials

if !isdir(Pkg.dir("TypedPolynomials"))
    println("Installing TypedPolynomials...")
    Pkg.clone("https://github.com/rdeits/TypedPolynomials.jl.git")
end
@reexport using TypedPolynomials

include("sphericalHarmonics.jl")
export ylm, rlm

include("sphericalHarmonicsExpansion.jl")
export sphericalHarmonicsExpansion, SphericalHarmonicCoefficients

include("fastfunc.jl")
export @fastfunc

include("sphericalQuadrature.jl")
export sphericalQuadrature

include("translation.jl")
export translation, translateRlm

end # module
