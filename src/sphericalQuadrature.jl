
function discreteIntegration(values::Array{Float64,1}, coordinates::Array{Float64,2}, l::Int, m::Int, x::PolyVar{true}, y::PolyVar{true}, z::PolyVar{true})

  sum = 0
  p = ylm(l,m,x,y,z)

  for k in 1:length(values)
    sum += values[k] * p([coordinates[k,1],coordinates[k,2],coordinates[k,3]], [x,y,z])
  end
  return sum * ((4*pi)/length(values))
end


function sphericalQuadrature(values::Array{Float64,1}, coordinates::Array{Float64,2}, L)

  @polyvar x y z
  C = [discreteIntegration(values,coordinates,l,m,x,y,z) for l=0:L for m=-l:l]
  return SphericalHarmonicCoefficients(C)
end
