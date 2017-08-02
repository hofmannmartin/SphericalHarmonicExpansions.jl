module SphericalHarmonics

using MultivariatePolynomials

#normalization factor
function sp_K(l::Int64, m::Int64)

  k::Float64 = 1
  for i in (l-m+1):(l+m)
    k *= i
  end

  return sqrt((2*l+1) / (4*pi*k))
end

#prefactor for the associated legendre polynomials
function sp_V(l::Int64, m::Int64)

  return ((-1)^m) / ((2^l)*factorial(l))
end

#trigonometric addition theorem for cos(m*phi)*sin(teta)^m
function sp_Tcos(x::MultivariatePolynomials.PolyVar{true}, y::MultivariatePolynomials.PolyVar{true}, m::Int64)

  sum = 0
  for j::Int64 in 0:floor(m/2)
    sum += ((-1)^j)*binomial(m, 2*j)*(y^(2*j))*(x^(m-2*j))
  end
  return sum
end

#trigonometric addition theorem for sin(m*phi)*sin(teta)^m
function sp_Tsin(x::MultivariatePolynomials.PolyVar{true}, y::MultivariatePolynomials.PolyVar{true}, m::Int64)

  sum = 0
  for j::Int64 in 0:floor((m-1)/2)
    sum += ((-1)^j)*binomial(m, 2*j + 1)*(y^(2*j + 1))*(x^(m-2*j-1))
  end
  return sum
end

#max sp_Y(20,20). For greater l there is an overflow error in sp_V(l,m)
#Tested up to l=4
#There is a factor (-1)^m, different from the table on wikipedia
function sp_Y(l::Int64, m::Int64, x::MultivariatePolynomials.PolyVar{true}, y::MultivariatePolynomials.PolyVar{true}, z::MultivariatePolynomials.PolyVar{true})

  if abs(m) > l
    println("Absolute Value of m can not be larger than l!")
  else

    p = (z^2 - 1)^l

    for i = 1:l+abs(m)
      p = differentiate(p, z)
    end

    if m > 0
      return sqrt(2)*sp_K(l, m)*sp_Tcos(x, y, m)*sp_V(l, m)*p

    elseif m < 0
      return sqrt(2)*sp_K(l, abs(m))*sp_Tsin(x, y, abs(m))*sp_V(l, abs(m))*p

    elseif m == 0
      return sp_K(l, 0)*sp_V(l, 0)*p

    else
      println("Error!")
    end
  end
end

#Returns a Vector with all spherical harmonics of l
function sp_Y(l::Int64, x::MultivariatePolynomials.PolyVar{true}, y::MultivariatePolynomials.PolyVar{true}, z::MultivariatePolynomials.PolyVar{true})

  Yl = [sp_Y(l,m,x,y,z) for m in -l:l]
  return Yl
end

#Returns all spherical harmonics up to l
function sp_allY(l::Int64, x::MultivariatePolynomials.PolyVar{true}, y::MultivariatePolynomials.PolyVar{true}, z::MultivariatePolynomials.PolyVar{true})

  Y = [sp_Y(ln,x,y,z) for ln in 0:l]
  return Y
end

end # module
