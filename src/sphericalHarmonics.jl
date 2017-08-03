#normalization factor
function ylmKCoefficient(l::Int64, m::Int64)

  k::Float64 = 1
  for i in (l-m+1):(l+m)
    k *= i
  end

  return sqrt((2*l+1) / (4*pi*k))
end

#
function ylmCosSinPolynomial(m::Int64, x::PolyVar{true}, y::PolyVar{true})

  sum = 0
  for j::Int64 in 0:floor(m/2)
    sum += ((-1)^j)*binomial(m, 2*j)*(y^(2*j))*(x^(m-2*j))
  end
  return sum
end

function ylmSinSinPolynomial(m::Int64, x::PolyVar{true}, y::PolyVar{true})

  sum = 0
  for j::Int64 in 0:floor((m-1)/2)
    sum += ((-1)^j)*binomial(m, 2*j + 1)*(y^(2*j + 1))*(x^(m-2*j-1))
  end
  return sum
end


function ylm(l::Int64, m::Int64, x::PolyVar{true}, y::PolyVar{true}, z::PolyVar{true})

  if abs(m) > l
    println("-l <= m <= l expected, but m = $m and l = $l.")
    error(BoundsError)
  else

    p = (z^2 - 1)^l

    for i = 1:l+abs(m)
      c = i <= l ? 1/(2*i) : -1.0
      p = c*differentiate(p, z)
    end

    if m > 0
      return sqrt(2)*ylmKCoefficient(l, m)*ylmCosSinPolynomial(m,x,y)*p

    elseif m < 0
      return sqrt(2)*ylmKCoefficient(l, abs(m))*ylmSinSinPolynomial(abs(m),x,y)*p

    elseif m == 0
      return ylmKCoefficient(l, 0)*p

    else
      println("Error!")
    end
  end
end
