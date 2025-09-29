"""
    rotation(c::SphericalHarmonicCoefficients, rot::Array{Float64,1})

*Input:*
- `c::SphericalHarmonicCoefficients`: Spherical harmonic coefficients
- `rot::Array{Float64,1}`: Euler angles for the rotation (ZYZ convention)

*Output:*
- `c_rot::SphericalHarmonicCoefficients`: Rotated spherical harmonic coefficients
"""
function rotation(shc::SphericalHarmonicCoefficients, rot::Array{Float64,1})
  shc_rot = deepcopy(shc)

  # run through all coefficients
  for l in 0:shc.L
    # calculate the new coefficients
    for m in -l:l
      shc_rot[l,m] = sum([shc[l,μ] * getRotation(l,m,μ,rot...) for μ = -l:l])
    end
  end

  return shc_rot
end


"""
    getRotation(l::Int,μ::Int,m::Int,α::Real,β::Real,γ::Real)

Get the elements `R^l_{μ,m}` of the rotation matrix for real spherical harmonics by Collado et al. (1989) (eqs. (10) - (18)), where (α,β,γ) are the Euler angles in the ZYZ convention.
"""
function getRotation(l::Int,μ::Int,m::Int,α::Real,β::Real,γ::Real)
  if μ == 0
    if m == 0
      # (14)
      return epsWignerdjmn(l,0,0,β)
    elseif m > 0
      # (13)
      return sqrt(2)*epsWignerdjmn(l,m,0,β) * cos(m*γ)
    else # m < 0
      # (15)
      return -sqrt(2)*epsWignerdjmn(l,-m,0,β) * sin(-m*γ)
    end
  elseif μ > 0
    if m == 0
      # (11)
      return sqrt(2)*epsWignerdjmn(l,0,μ,β) * cos(μ*α)
    elseif m > 0
      # (10)
      return epsWignerdjmn(l,-μ,-m,β) * cos(μ*α + m*γ) + (-1)^μ * epsWignerdjmn(l,μ,-m,β) * cos(μ*α - m*γ)
    else # m < 0
      # (12)
      return -epsWignerdjmn(l,-μ,m,β) * sin(μ*α - m*γ) + (-1)^μ * epsWignerdjmn(l,μ,m,β) * sin(μ*α + m*γ)
    end
  else # μ < 0
    if m == 0
      # (17)
      return sqrt(2)*epsWignerdjmn(l,0,-μ,β) * sin(-μ*α)
    elseif m > 0
      # (16)
      return epsWignerdjmn(l,μ,-m,β) * sin(-μ*α + m*γ) + (-1)^μ * epsWignerdjmn(l,-μ,-m,β) * sin(-μ*α - m*γ)
    else # m < 0
      # (18)
      return epsWignerdjmn(l,μ,m,β) * cos(-μ*α - m*γ) - (-1)^μ * epsWignerdjmn(l,-μ,m,β) * cos(-μ*α + m*γ)
    end
  end
end

"""
Set values smaller than the machine epsilon to zero.
"""
function epsWignerdjmn(l,μ,m,β)
  value = WignerD.wignerdjmn(l,μ,m,β)
  return abs(value) < eps() ? 0.0 : value
end


"""
    pointReflection(c::SphericalHarmonicCoefficients)

Apply a point reflection to the spherical harmonic coefficients at the origin.

*Input:*
- `c::SphericalHarmonicCoefficients`: Spherical harmonic coefficients

*Output:*
- `c_mir::SphericalHarmonicCoefficients`: Reflected spherical harmonic coefficients
"""
function pointReflection(shc::SphericalHarmonicCoefficients)
  shc_pr = deepcopy(shc)
    
  # run through all coefficients and calculate the new coefficients
  for l in 0:shc.L
    for m in -l:l
      shc_pr[l,m] = (-1)^l * shc[l,m]
    end
  end

  return shc_pr
end