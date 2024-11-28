#normalization factor
function ylmCoefficient(l, m)
  kpi = 4*pi
  for i in l-m+1:l+m
      kpi *= i
  end
  out = sqrt((2*l+1)/kpi)
  if out!=0
    return sqrt((2*l+1)/kpi)
  else
    error("ylmCoefficient could not be represented in floating point precicion")
  end    
end

function ylmCosSinPolynomial(m, x, y)
	terms = [((-1)^j)*binomial(m, 2*j)*y^(2*j)*x^(m-2*j) for j in 0:floor(Int,m/2)]
	return polynomial(terms)
end

function ylmSinSinPolynomial(m, x, y)
	terms = [((-1)^j)*binomial(m, 2*j + 1)*y^(2*j + 1)*x^(m-2*j-1) for j in 0:floor(Int,(m-1)/2)]
	return polynomial(terms)
end

"""
    ylm(l::Int64, m::Int64, x::Variable, y::Variable, z::Variable)
*Description:*  Calculation of the spherical harmonic for a given order (l,m) in Cartesian coordinates\\

*Input:*  `l`       - Degree of the spherical harmonic\\
          `m`       - Order of the spherical harmonic\\
          `x, y, z` - Cartesian coordinates\\

*Output:*  Spherical harmonic polynomial
"""
function ylm(l, m, x, y, z)
	if abs(m) > l
		throw(DomainError(m,"-l <= m <= l expected, but m = $m and l = $l."))
	end

	terms = [Float64((-1)^k*binomial(l,k))*z^(2*(l-k)) for k=0:ceil(Int,(l-abs(m))/2)]
	p = polynomial(terms) + 0.0*(x+y)
	for i = 1:l+abs(m)
		c = i <= l ? 1/(2*i) : 1.0
		p = c*differentiate(p, z, Val{1}())
	end

	if m > 0
		return sqrt(2)*ylmCoefficient(l, m)*ylmCosSinPolynomial(m,x,y)*p
	elseif m < 0
		return sqrt(2)*ylmCoefficient(l, abs(m))*ylmSinSinPolynomial(abs(m),x,y)*p
	else
		return ylmCoefficient(l, 0)*p
	end
end

# expand (x²+y²+z²)^n
function trinomialExpansion(n , x, y, z)
	multiindices = [(i,j,n-i-j) for i in 0:n for j in 0:n-i]
	terms = [multinomial(i,j,k)*x^(2*i)*y^(2*j)*z^(2*k) for (i,j,k) in multiindices]
	return polynomial(terms)
end

# multiplying r^l*ylm(x,y,z)
function rlylm(l, m, x, y, z)
	p = ylm(l,m,x,y,z)

  pout = 0.0*(x+y+z)
	for t in terms(p)
		deg = degree(t) # Gibt den gesamten Grad des Monoms an
		degR = l-deg # durch das Kürzen ergibt sich ein Grad von l-deg fuer r
		pout += trinomialExpansion(div(degR,2),x,y,z)*t # r² wird durch x²+y²+z² ersetzt
	end

	return pout
end

# solid harmonics
function zlm(l, m, x, y, z)
	rlm = rlylm(l,m,x,y,z)
	rlm = sqrt(4*pi/(2*l+1))*rlm
	return rlm
end
