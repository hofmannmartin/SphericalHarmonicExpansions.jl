
"""
    sphericalHarmonicsExpansion(Clm::Array{Float64,1}, x::PolyVar{true}, y::PolyVar{true}, z::PolyVar{true})
*Description:*  Calculation of the spherical harmonics expansion in Cartesian coordinates
                for given coefficients which define the maximum degree of the spherical harmonics\\
*Input:*  `Clm`       - Array with coefficients (length = (l+1)², l = max. deg. of the spherical harmonics)\\
          `x, y, z`   - Cartesian coordinates\\
*Output:*  Spherical harmonics expansion
"""
function sphericalHarmonicsExpansion(Clm::Array{Float64,1}, x::PolyVar{true}, y::PolyVar{true}, z::PolyVar{true})

  N = size(Clm,1)
  L = isqrt(N-1)
  realN = L^2+2L+1
  if N != realN
    println("N = l²+2l+1 expected, N = $N but should be N = $realN for l = $L !")
    error(BoundsError)
  else

    sum = 0

    for l in 0:L
      for m in -l:l
        if Clm[l*(l+1)+m+1] != 0
          sum += Clm[l*(l+1)+m+1] * ylm(l,m,x,y,z)
        end
      end
    end
    return sum
  end
end
