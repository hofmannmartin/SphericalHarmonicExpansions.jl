#normalization factor
function ylmKCoefficient(l::Int64, m::Int64)

  k::Float64 = 1
  for i in (l-m+1):(l+m)
    k *= i
  end

  return sqrt((2*l+1) / (4*pi*k))
end

function ylmCosSinPolynomial(m::Int64, x::Variable, y::Variable)

  sum = 0.0*(x+y)
  for j::Int64 in 0:floor(m/2)
    sum += ((-1)^j)*binomial(m, 2*j)*(y^(2*j))*(x^(m-2*j))
  end
  return sum
end

function ylmSinSinPolynomial(m::Int64, x::Variable, y::Variable)

  sum = 0.0*(x+y)
  for j::Int64 in 0:floor((m-1)/2)
    sum += ((-1)^j)*binomial(m, 2*j + 1)*(y^(2*j + 1))*(x^(m-2*j-1))
  end
  return sum
end

"""
    ylm(l::Int64, m::Int64, x::Variable, y::Variable, z::Variable)
*Description:*  Calculation of the spherical harmonic for a given order (l,m) in Cartesian coordinates\\

*Input:*  `l`       - Degree of the spherical harmonic\\
          `m`       - Order of the spherical harmonic\\
          `x, y, z` - Cartesian coordinates\\

*Output:*  Spherical harmonic polynomial
"""
function ylm(l::Int64, m::Int64, x::Variable, y::Variable, z::Variable)

  if abs(m) > l
    throw(DomainError(m,"-l <= m <= l expected, but m = $m and l = $l."))
  end

  p = (z^2 - 1)^l + 0.0*(x+y)

  for i = 1:l+abs(m)
    c = i <= l ? 1/(2*i) : 1.0
    p = c*differentiate(p, z)
  end

  if m > 0
    return sqrt(2)*ylmKCoefficient(l, m)*ylmCosSinPolynomial(m,x,y)*p
  elseif m < 0
    return sqrt(2)*ylmKCoefficient(l, abs(m))*ylmSinSinPolynomial(abs(m),x,y)*p
  else
    return ylmKCoefficient(l, 0)*p
  end
end

# multiplying r^l*ylm(x,y,z)
function rlylm(l::Int64, m::Int, x::Variable, y::Variable, z::Variable)
	p = ylm(l,m,x,y,z)
	tout = []
	# Zerlegung des Polynoms in Terme:
	for t in terms(p)
		deg = degree(monomial(t)) # Gibt den gesamten Grad des Monoms an
		degR = l-deg # durch das Kürzen ergibt sich ein Grad von l-deg fuer r
		push!(tout,(x^2+y^2+z^2)^Int(degR/2)*t) # r² wird durch x²+y²+z² ersetzt
	end

	return polynomial(tout)
end

# solid harmonics
function zlm(l::Int64, m::Int64, x::Variable, y::Variable, z::Variable)
	rlm = rlylm(l,m,x,y,z)
	rlm = sqrt(4*pi/(2*l+1))*rlm
	return rlm
end
