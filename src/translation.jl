"""
    translation(C::SphericalHarmonicCoefficients,v::Array{Float64,1})
*Description:* Translation of the coefficients: Shifting the expansion point by v\\
\\
*Input:*
- `C`  - Coefficients
- `v`  - shift vector (length(v) = 3)
*Output:*
- `cT` - Translated coefficients, type: SphericalHarmonicCoefficients (cT.R = C.R, cT.solid = C.solid)

"""
function translation(C::SphericalHarmonicCoefficients,v::Array{Float64,1})

    @polyvar x y z

    cT = deepcopy(C) # translated coefficients

    # turn spherical coefficients into solid coefficients
    if !(C.solid)
        solid!(C)
    end

    # [l,0] for l >= 0 (eq. (24) & (25))
    for l=0:C.L
      cT[l,0] = 0
      for λ=l:C.L
        for μ = -(λ-l):λ-l
          p = zlm(λ-l,μ,x,y,z)
          cT[l,0] += C[λ,μ] * σ(λ,μ,l,0,1) * p((x,y,z) => v)
        end
      end
    end

    # [l,m] for l ≠ 0, m > 0 (eq. (22))
    for l=1:C.L
      for m=1:l
        cT[l,m] = 0
        # 1: Summand from (13)
        for λ=l:C.L
          for μ = m:m-(l-λ)
            p = zlm(λ-l,μ-m,x,y,z)
            cT[l,m] += C[λ,μ] * σ(λ,μ,l,m,1) * p((x,y,z) => v)
          end
        end
        # 2: Summand from (13)
        for λ=l+1:C.L
          for μ = max(1,m-(λ-l)):m-1
            p = zlm(λ-l,m-μ,x,y,z)
            cT[l,m] += C[λ,μ] * σ(λ,μ,l,m,2) * p((x,y,z) => v)
          end
        end
        # 3: Summand from (13)
        for λ=l+1:C.L
          for μ = 1:-m-(l-λ)
            p = zlm(λ-l,μ+m,x,y,z)
            cT[l,m] += C[λ,μ] * σ(λ,μ,l,-m,3) * p((x,y,z) => v)
          end
        end
        # 4: Summand from (14)
        for λ=m+l:C.L
          p = zlm(λ-l,m,x,y,z)
          cT[l,m] += C[λ,0] * (σ(λ,0,l,m,2) + σ(λ,0,l,-m,3)) * p((x,y,z) => v)
        end
        # 5: Summand from (15)
        for λ=l:C.L
          for μ = -m-(λ-l):-m-1 #min(-1,-m)
            p = zlm(λ-l,μ+m,x,y,z)
            cT[l,m] += C[λ,μ] * σ(λ,μ,l,-m,1) * p((x,y,z) => v)
          end
        end
        # 6: Summand from (15)
        for λ=l+1:C.L
          for μ = -m+1:min(-1,-m-(l-λ))
            p = zlm(λ-l,-(μ+m),x,y,z)
            cT[l,m] -= C[λ,μ] * σ(λ,μ,l,-m,2) * p((x,y,z) => v)
          end
        end
        # 7: Summand from (15)
        for λ=l+1:C.L
          for μ = m-(λ-l):-1
            p = zlm(λ-l,μ-m,x,y,z)
            cT[l,m] += C[λ,μ] * σ(λ,μ,l,m,3) * p((x,y,z) => v)
          end
        end
      end
    end

    # [l,m] for l ≠ 0, m < 0 (eq. (23))
    for l=1:C.L
      for m=-l:-1
        cT[l,m] = 0
        # 1: Summand from (13)
        for λ=l:C.L
          for μ = -m+1:-m-(l-λ)
            p = zlm(λ-l,-(μ+m),x,y,z)
            cT[l,m] -= C[λ,μ] * σ(λ,μ,l,-m,1) * p((x,y,z) => v)
          end
        end
        # 2: Summand from (13)
        for λ=l+1:C.L
          for μ = max(1,-m-(λ-l)):-m-1
            p = zlm(λ-l,μ+m,x,y,z)
            cT[l,m] += C[λ,μ] * σ(λ,μ,l,-m,2) * p((x,y,z) => v)
          end
        end
        # 3: Summand from (13)
        for λ=l+1:C.L
          for μ = 1:m-(l-λ)
            p = zlm(λ-l,-(μ-m),x,y,z)
            cT[l,m] += C[λ,μ] * σ(λ,μ,l,m,3) * p((x,y,z) => v)
          end
        end
        # 4: Summand from (14)
        for λ=l-m:C.L
          p = zlm(λ-l,m,x,y,z)
          cT[l,m] += C[λ,0] * (σ(λ,0,l,-m,2) + σ(λ,0,l,m,3)) * p((x,y,z) => v)
        end
        # 5: Summand from (15)
        for λ=l:C.L
          for μ = m-(λ-l):m#-1 #min(-1,m)
            p = zlm(λ-l,m-μ,x,y,z)
            cT[l,m] += C[λ,μ] * σ(λ,μ,l,m,1) * p((x,y,z) => v)
          end
        end
        # 6: Summand from (15)
        for λ=l+1:C.L
          for μ = m+1:min(-1,m-(l-λ))
            p = zlm(λ-l,μ-m,x,y,z)
            cT[l,m] += C[λ,μ] * σ(λ,μ,l,m,2) * p((x,y,z) => v)
          end
        end
        # 7: Summand from (15)
        for λ=l+1:C.L
          for μ = -m-(λ-l):-1
            p = zlm(λ-l,-(μ+m),x,y,z)
            cT[l,m] -= C[λ,μ] * σ(λ,μ,l,-m,3) * p((x,y,z) => v)
          end
        end
      end
    end

    # turn solid coefficients into spherical coefficients
    if !(cT.solid)
        cT.solid = true
        spherical!(cT)
        spherical!(C)
    end

    return cT
end

## prefactors
function σ(l,m,λ,μ,num::Int)
  σ = sqrt((factorial(l+m) * factorial(l-m)) / (factorial(λ+μ) * factorial(λ-μ) * factorial(l-λ+m-μ) * factorial(l-λ-m+μ)))
  if num == 1
    σ *= (μ != 0 && μ != m && m != 0) ? 1/sqrt(2) : 1
  else
    # num = 2
    σ *= (-1)^(μ-m) / 2
    σ *= (m != 0) ? sqrt(2) : 1
    # num = 3
    σ *= (num == 3) ? (-1)^m : 1
  end
  return σ
end

function translateRlm(l::Int64, m::Int64,vx,vy,vz)

    @polyvar x y z

    Rlmt = 0;

    # two cases m≥0 and m<0, each with three summands:
    if m >= 0
    for λ=0:l
        for μ=max(0,λ-l+m):min(λ,m)
            if m == 0
                p = zlm(l-λ,0,x,y,z);
                mult = 1/(factorial(λ)*factorial(l-λ));
                sum1 = zlm(λ,0,x,y,z)*p((x,y,z)=>(vx,vy,vz));
                Rlmt += mult*sum1;
            elseif μ == 0 # => m-μ = m > 0
                p = zlm(l-λ,m,x,y,z);
                mult = 1/(factorial(λ)*sqrt(2*factorial(l-λ+m)*factorial(l-λ-m)));
                sum1 = zlm(λ,0,x,y,z)*p((x,y,z)=>(vx,vy,vz));
                Rlmt += mult*sum1;
            elseif μ == m # => m-μ = 0
                p = zlm(l-λ,0,x,y,z);
                mult = 1/(factorial(l-λ)*sqrt(2*factorial(λ+m)factorial(λ-m)));
                sum1 = zlm(λ,m,x,y,z)*p((x,y,z)=>(vx,vy,vz));
                Rlmt += mult*sum1;
            else # m-μ > 0, μ > 0
                p = zlm(l-λ,m-μ,x,y,z);
                q = zlm(l-λ,-(m-μ),x,y,z);
                mult = 1/sqrt(4*factorial(λ+μ)*factorial(λ-μ)*factorial(l-λ+m-μ)*factorial(l-λ-(m-μ)));
                sum1 = zlm(λ,μ,x,y,z)*p((x,y,z)=>(vx,vy,vz));
                sum2 = zlm(λ,-μ,x,y,z)*q((x,y,z)=>(vx,vy,vz));
                Rlmt += mult*(sum1 - sum2);
            end
        end
    end

    for λ=m+1:l-1
        for μ=m+1:min(λ,-λ+l+m) # μ > 0, m-μ < 0
            p = zlm(l-λ,abs(m-μ),x,y,z);
            q = zlm(l-λ,-abs(m-μ),x,y,z);
            mult = (-1)^(m-μ)/sqrt(4*factorial(λ+μ)*factorial(λ-μ)*factorial(l-λ+(m-μ))*factorial(l-λ-(m-μ)));
            sum1 = zlm(λ,μ,x,y,z)*p((x,y,z)=>(vx,vy,vz));
            sum2 = zlm(λ,-μ,x,y,z)*q((x,y,z)=>(vx,vy,vz));
            Rlmt += mult*(sum1 + sum2);
        end
    end

    for λ=1:l-m-1
        for μ=max(-λ,λ-l+m):-1 # μ < 0, m-μ > 0
            p = zlm(l-λ,m-μ,x,y,z);
            q = zlm(l-λ,-(m-μ),x,y,z);
            mult = (-1)^(μ)/sqrt(4*factorial(λ+μ)*factorial(λ-μ)*factorial(l-λ+m-μ)*factorial(l-λ-(m-μ)));
            sum1 = zlm(λ,abs(μ),x,y,z)*p((x,y,z)=>(vx,vy,vz));
            sum2 = zlm(λ,-abs(μ),x,y,z)*q((x,y,z)=>(vx,vy,vz));
            Rlmt += mult*(sum1 + sum2);
        end
    end

    Rlmt *= sqrt(factorial(l-m)*factorial(l+m));
    Rlmt *= (m == 0) ? 1 : sqrt(2);

    else # m < 0
        for λ=0:l
            for μ=max(m,-λ):min(0,-λ+l+m)
                if μ == 0 # => m-μ = m < 0
                    p = zlm(l-λ,-abs(m),x,y,z);
                    mult = (-1)^m/(factorial(λ)*sqrt(2*factorial(l-λ+m)*factorial(l-λ-m)));
                    sum1 = zlm(λ,0,x,y,z)*p((x,y,z)=>(vx,vy,vz));
                    Rlmt += (-1)*mult*sum1;
                elseif μ == m # => m-μ = 0
                    p = zlm(l-λ,0,x,y,z);
                    mult = (-1)^m/(factorial(l-λ)*sqrt(2*factorial(λ+m)factorial(λ-m)));
                    sum1 = zlm(λ,-abs(m),x,y,z)*p((x,y,z)=>(vx,vy,vz));
                    Rlmt += (-1)*mult*sum1;
                else # μ < 0, m-μ < 0
                    p = zlm(l-λ,abs(m-μ),x,y,z);
                    q = zlm(l-λ,-abs(m-μ),x,y,z);
                    mult = (-1)^m/sqrt(4*factorial(λ+μ)*factorial(λ-μ)*factorial(l-λ+(m-μ))*factorial(l-λ-(m-μ)));
                    sum1 = zlm(λ,-abs(μ),x,y,z)*p((x,y,z)=>(vx,vy,vz));
                    sum2 = zlm(λ,abs(μ),x,y,z)*q((x,y,z)=>(vx,vy,vz));
                    Rlmt += (-1)*mult*(sum1 + sum2);
                end
            end
        end

        for λ=-m+1:l-1
            for μ=max(-λ,λ-l+m):m-1 # μ < 0, m-μ > 0
                p = zlm(l-λ,-(m-μ),x,y,z);
                q = zlm(l-λ,(m-μ),x,y,z);
                mult = (-1)^(μ)/sqrt(4*factorial(λ+μ)*factorial(λ-μ)*factorial(l-λ+m-μ)*factorial(l-λ-(m-μ)));
                sum1 = zlm(λ,abs(μ),x,y,z)*p((x,y,z)=>(vx,vy,vz));
                sum2 = zlm(λ,-abs(μ),x,y,z)*q((x,y,z)=>(vx,vy,vz));
                Rlmt += mult*(sum1 - sum2);
            end
        end

        for λ=1:l+m-1
            for μ=1:min(λ,-λ+l+m) # μ > 0, m-μ < 0
                p = zlm(l-λ,abs(m-μ),x,y,z);
                q = zlm(l-λ,-abs(m-μ),x,y,z);
                mult = (-1)^(m-μ)/sqrt(4*factorial(λ+μ)*factorial(λ-μ)*factorial(l-λ+(m-μ))*factorial(l-λ-(m-μ)));
                sum1 = zlm(λ,-μ,x,y,z)*p((x,y,z)=>(vx,vy,vz));
                sum2 = zlm(λ,μ,x,y,z)*q((x,y,z)=>(vx,vy,vz));
                Rlmt += mult*(sum1 - sum2);
            end
        end

        Rlmt *= (-1)*(-1)^m*sqrt(2*factorial(l-abs(m))*factorial(l+abs(m)));
    end
    return Rlmt
end

# Error propagation
function errorTranslation(C::SphericalHarmonicCoefficients,v)

    @polyvar x y z

    cT = deepcopy(C) # translated coefficients

    # turn spherical coefficients into solid coefficients
    if !(C.solid)
        solid!(C)
    end

    # [l,0] for l >= 0 (eq. (24) & (25))
    for l=0:C.L
      cT[l,0] = 0
      for λ=l:C.L
        for μ = -(λ-l):λ-l
          p = zlm(λ-l,μ,x,y,z)
          cT[l,0] += (C[λ,μ] * σ(λ,μ,l,0,1) * p((x,y,z) => v))^2
        end
      end
    end

    # [l,m] for l ≠ 0, m > 0 (eq. (22))
    for l=1:C.L
      for m=1:l
        cT[l,m] = 0
        # 1: Summand from (13)
        for λ=l:C.L
          for μ = m:m-(l-λ)
            p = zlm(λ-l,μ-m,x,y,z)
            cT[l,m] += (C[λ,μ] * σ(λ,μ,l,m,1) * p((x,y,z) => v))^2
          end
        end
        # 2: Summand from (13)
        for λ=l+1:C.L
          for μ = max(1,m-(λ-l)):m-1
            p = zlm(λ-l,m-μ,x,y,z)
            cT[l,m] += (C[λ,μ] * σ(λ,μ,l,m,2) * p((x,y,z) => v))^2
          end
        end
        # 3: Summand from (13)
        for λ=l+1:C.L
          for μ = 1:-m-(l-λ)
            p = zlm(λ-l,μ+m,x,y,z)
            cT[l,m] += (C[λ,μ] * σ(λ,μ,l,-m,3) * p((x,y,z) => v))^2
          end
        end
        # 4: Summand from (14)
        for λ=m+l:C.L
          p = zlm(λ-l,m,x,y,z)
          cT[l,m] += (C[λ,0] * (σ(λ,0,l,m,2) + σ(λ,0,l,-m,3)) * p((x,y,z) => v))^2
        end
        # 5: Summand from (15)
        for λ=l:C.L
          for μ = -m-(λ-l):-m-1 #min(-1,-m)
            p = zlm(λ-l,μ+m,x,y,z)
            cT[l,m] += (C[λ,μ] * σ(λ,μ,l,-m,1) * p((x,y,z) => v))^2
          end
        end
        # 6: Summand from (15)
        for λ=l+1:C.L
          for μ = -m+1:min(-1,-m-(l-λ))
            p = zlm(λ-l,-(μ+m),x,y,z)
            cT[l,m] -= (C[λ,μ] * σ(λ,μ,l,-m,2) * p((x,y,z) => v))^2
          end
        end
        # 7: Summand from (15)
        for λ=l+1:C.L
          for μ = m-(λ-l):-1
            p = zlm(λ-l,μ-m,x,y,z)
            cT[l,m] += (C[λ,μ] * σ(λ,μ,l,m,3) * p((x,y,z) => v))^2
          end
        end
      end
    end

    # [l,m] for l ≠ 0, m < 0 (eq. (23))
    for l=1:C.L
      for m=-l:-1
        cT[l,m] = 0
        # 1: Summand from (13)
        for λ=l:C.L
          for μ = -m+1:-m-(l-λ)
            p = zlm(λ-l,-(μ+m),x,y,z)
            cT[l,m] -= (C[λ,μ] * σ(λ,μ,l,-m,1) * p((x,y,z) => v))^2
          end
        end
        # 2: Summand from (13)
        for λ=l+1:C.L
          for μ = max(1,-m-(λ-l)):-m-1
            p = zlm(λ-l,μ+m,x,y,z)
            cT[l,m] += (C[λ,μ] * σ(λ,μ,l,-m,2) * p((x,y,z) => v))^2
          end
        end
        # 3: Summand from (13)
        for λ=l+1:C.L
          for μ = 1:m-(l-λ)
            p = zlm(λ-l,-(μ-m),x,y,z)
            cT[l,m] += (C[λ,μ] * σ(λ,μ,l,m,3) * p((x,y,z) => v))^2
          end
        end
        # 4: Summand from (14)
        for λ=l-m:C.L
          p = zlm(λ-l,m,x,y,z)
          cT[l,m] += (C[λ,0] * (σ(λ,0,l,-m,2) + σ(λ,0,l,m,3)) * p((x,y,z) => v))^2
        end
        # 5: Summand from (15)
        for λ=l:C.L
          for μ = m-(λ-l):m#-1 #min(-1,m)
            p = zlm(λ-l,m-μ,x,y,z)
            cT[l,m] += (C[λ,μ] * σ(λ,μ,l,m,1) * p((x,y,z) => v))^2
          end
        end
        # 6: Summand from (15)
        for λ=l+1:C.L
          for μ = m+1:min(-1,m-(l-λ))
            p = zlm(λ-l,μ-m,x,y,z)
            cT[l,m] += (C[λ,μ] * σ(λ,μ,l,m,2) * p((x,y,z) => v))^2
          end
        end
        # 7: Summand from (15)
        for λ=l+1:C.L
          for μ = -m-(λ-l):-1
            p = zlm(λ-l,-(μ+m),x,y,z)
            cT[l,m] -= (C[λ,μ] * σ(λ,μ,l,-m,3) * p((x,y,z) => v))^2
          end
        end
      end
    end

    cT.c = sqrt.(cT.c)

    # turn solid coefficients into spherical coefficients
    if !(cT.solid)
        cT.solid = true
        spherical!(cT)
        spherical!(C)
    end

    return cT
end
