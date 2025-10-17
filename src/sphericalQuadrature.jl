
function discreteIntegration(values::AbstractArray{<:Real,1}, coordinates::AbstractArray{<:Real,2}, l::Int, m::Int)

  @polyvar x y z

  sum = 0
  p = ylm(l,m,x,y,z)

  for k in eachindex(values)
    sum += values[k] * p((x,y,z)=>(coordinates[k,1],coordinates[k,2],coordinates[k,3]))
  end
  return sum * ((4*pi)/length(values))
end


function sphericalQuadrature(values::AbstractArray{<:Real,1}, coordinates::AbstractArray{<:Real,2}, L)

  C = [discreteIntegration(values,coordinates,l,m) for l=0:L for m=-l:l]
  return SphericalHarmonicCoefficients(C)
end

# Error propagation
function errorDiscreteIntegration(values::AbstractArray{<:Real,1}, coordinates::AbstractArray{<:Real,2}, l::Int, m::Int)

  @polyvar x y z

  sum = 0
  p = ylm(l,m,x,y,z)

  for k in eachindex(values)
    sum += (values[k] * p((x,y,z)=>(coordinates[k,1],coordinates[k,2],coordinates[k,3])))^2
  end
  return sqrt(sum) * ((4*pi)/length(values))
end


function errorSphericalQuadrature(values::AbstractArray{<:Real,1}, coordinates::AbstractArray{<:Real,2}, L)

  C = [errorDiscreteIntegration(values,coordinates,l,m) for l=0:L for m=-l:l]
  return SphericalHarmonicCoefficients(C)
end
