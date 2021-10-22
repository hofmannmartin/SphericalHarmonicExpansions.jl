
function discreteIntegration(values::Array{Float64,1}, coordinates::AbstractArray{Float64,2}, l::Int, m::Int)

  @polyvar x y z

  sum = 0
  p = ylm(l,m,x,y,z)

  for k in 1:length(values)
    sum += values[k] * p((x,y,z)=>(coordinates[k,1],coordinates[k,2],coordinates[k,3]))
  end
  return sum * ((4*pi)/length(values))
end


function sphericalQuadrature(values::Array{Float64,1}, coordinates::AbstractArray{Float64,2}, L)

  C = [discreteIntegration(values,coordinates,l,m) for l=0:L for m=-l:l]
  return SphericalHarmonicCoefficients(C)
end

# Error propagation
function errorDiscreteIntegration(values::Array{Float64,1}, coordinates::AbstractArray{Float64,2}, l::Int, m::Int)

  @polyvar x y z

  sum = 0
  p = ylm(l,m,x,y,z)

  for k in 1:length(values)
    sum += (values[k] * p((x,y,z)=>(coordinates[k,1],coordinates[k,2],coordinates[k,3])))^2
  end
  return sqrt(sum) * ((4*pi)/length(values)) * sqrt(4*pi/(2*l+1))
end


function errorSphericalQuadrature(values::Array{Float64,1}, coordinates::AbstractArray{Float64,2}, L)

  C = [errorDiscreteIntegration(values,coordinates,l,m) for l=0:L for m=-l:l]
  return SphericalHarmonicCoefficients(C)
end
