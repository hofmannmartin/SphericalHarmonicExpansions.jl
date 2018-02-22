import Base.setindex!, Base.getindex, Base.isapprox, Base.==, Base.!=,
        Base.+, Base.-, Base.*, Base./, Base.normalize

type SphericalHarmonicCoefficients
  c::Vector{T} where T<:Real
  L::Int

  function SphericalHarmonicCoefficients(c::Vector{T} where T<:Real)
    if (!isinteger(sqrt(length(c))) || length==0)
       println("input vector needs o be of size (L+1)², where L ∈ ℕ₀.")
       throw(DomainError())
    end

    L = convert(Int,sqrt(length(c)))-1
    return new(c,L)
  end

  function SphericalHarmonicCoefficients(L::Int)
    if L<0
      println("input vector needs o be of size (L+1)², where L ∈ ℕ₀.")
      throw(DomainError())
    end
    return new(zeros((L+1)^2),UInt(L))
  end
end

isapprox(shc1::SphericalHarmonicCoefficients, shc2::SphericalHarmonicCoefficients; kargs...) = isapprox(shc1.c,shc2.c;kargs...)
getindex(shc::SphericalHarmonicCoefficients,I) = getindex(shc.c,I)
getindex(shc::SphericalHarmonicCoefficients,l,m) = getindex(shc.c, l*(l+1)+m+1)
setindex!(shc::SphericalHarmonicCoefficients,value,I) = setindex!(shc.c,value,I)
setindex!(shc::SphericalHarmonicCoefficients,value,l,m) = setindex!(shc.c,value,l*(l+1)+m+1)
==(shca::SphericalHarmonicCoefficients, shcb::SphericalHarmonicCoefficients) = shca.c == shcb.c
!=(shca::SphericalHarmonicCoefficients, shcb::SphericalHarmonicCoefficients) = shca.c != shcb.c
+(shca::SphericalHarmonicCoefficients, shcb::SphericalHarmonicCoefficients) = SphericalHarmonicCoefficients(shca.c + shcb.c)
-(shca::SphericalHarmonicCoefficients, shcb::SphericalHarmonicCoefficients) = SphericalHarmonicCoefficients(shca.c - shcb.c)
*(shca::SphericalHarmonicCoefficients, shcb::SphericalHarmonicCoefficients) = SphericalHarmonicCoefficients(shca.c * shcb.c)
/(shca::SphericalHarmonicCoefficients, shcb::SphericalHarmonicCoefficients) = SphericalHarmonicCoefficients(shca.c / shcb.c)
+(shc::SphericalHarmonicCoefficients, value) = SphericalHarmonicCoefficients(shc.c .+ value)
-(shc::SphericalHarmonicCoefficients, value) = SphericalHarmonicCoefficients(shc.c .- value)
*(shc::SphericalHarmonicCoefficients, value) = SphericalHarmonicCoefficients(shc.c .* value)
/(shc::SphericalHarmonicCoefficients, value) = SphericalHarmonicCoefficients(shc.c ./ value)
normalize(shc::SphericalHarmonicCoefficients,R::Float64) = SphericalHarmonicCoefficients([shc[l,m] *= 1/(R^l) for l = 0:shc.L for m = -l:l])

"""
    sphericalHarmonicsExpansion(Clm::Array{Float64,1}, x::Variable, y::Variable, z::Variable)
*Description:*  Calculation of the spherical harmonics expansion in Cartesian coordinates
                for given coefficients which define the maximum degree of the spherical harmonics\\
*Input:*  `Clm`       - Array with coefficients (length = (l+1)², l = max. deg. of the spherical harmonics)\\
          `x, y, z`   - Cartesian coordinates\\
*Output:*  Spherical harmonics expansion
"""
function sphericalHarmonicsExpansion(Clm::SphericalHarmonicCoefficients, x::Variable, y::Variable, z::Variable)

  sum = 0

  for l in 0:Clm.L
    for m in -l:l
      if Clm[l,m] != 0
        sum += Clm[l,m] * rlylm(l,m,x,y,z)
      end
    end
  end
  return sum
end

"""
    solidHarmonicsExpansion(Clm::Array{Float64,1}, x::Variable, y::Variable, z::Variable)
*Description:*  Calculation of the solid harmonics expansion in Cartesian coordinates
                for given coefficients which define the maximum degree of the spherical harmonics\\
*Input:*  `Clm`       - Array with coefficients (length = (l+1)², l = max. deg. of the spherical harmonics)\\
          `x, y, z`   - Cartesian coordinates\\
*Output:*  Solid harmonics expansion
"""
function solidHarmonicsExpansion(Clm::SphericalHarmonicCoefficients, x::Variable, y::Variable, z::Variable)

  sum = 0

  for l in 0:Clm.L
    for m in -l:l
      if Clm[l,m] != 0
        sum += Clm[l,m] * rlm(l,m,x,y,z)
      end
    end
  end
  return sum
end
