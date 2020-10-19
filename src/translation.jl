"""
    translation(C::SphericalHarmonicCoefficients,v::Array{Float64,1})
*Description:* Translation of the coefficients: Shifting the expansion point by v\\
\\
*Input:*
- `C`  - Coefficients
- `v`  - shift vector (length(v) = 3)
*Output:*
- `cShifted` - Shifted coefficients, type: SphericalHarmonicCoefficients (cShifted.R = C.R, cShifted.solid = C.solid)

"""
function translation(C::SphericalHarmonicCoefficients,v::Array{Float64,1})

    @polyvar x y z

    cShifted = deepcopy(C)

    vx = v[1]
    vy = v[2]
    vz = v[3]

    # turn solid coefficients into spherical coefficients
    if C.solid
        spherical!(C)
    end

    # c[0,0]:
    sum = 0
    # Summand aus (I) - m > 0
    for λ=1:C.L
        for μ = 1:λ
            p = rlm(λ,μ,x,y,z)
            sum += ω₊(C[λ,μ],λ,μ,0,0,1)*p((x,y,z)=>(vx,vy,vz))
        end
    end
    # Summand aus (II) - m = 0
    for λ=0:C.L
        p = rlm(λ,0,x,y,z)
        sum += ω₀(C[λ,0],λ,0,0,1)*p((x,y,z)=>(vx,vy,vz))
    end
    # Summand aus (III) - m < 0
    for λ=1:C.L
        for μ = -λ:-1
            p = rlm(λ,μ,x,y,z)
            sum += ω₋(C[λ,μ],λ,μ,0,0,1)*p((x,y,z)=>(vx,vy,vz))
        end
    end
    cShifted[0,0] = sqrt(4*pi)*sum
    sum = 0

    # cShifted[l,0] (l > 0):
    for l=1:C.L
        # Summand aus (I) - m > 0
        for λ=l:C.L
            for μ = 1:λ-l
                p = rlm(λ-l,μ,x,y,z)
                sum += ω₊(C[λ,μ],λ,μ,l,0,1)*p((x,y,z)=>(vx,vy,vz))
            end
        end
        # Summand aus (II) - m = 0
        for λ=l:C.L
            p = rlm(λ-l,0,x,y,z)
            sum += ω₀(C[λ,0],λ,l,0,1)*p((x,y,z)=>(vx,vy,vz))
        end
        # Summand aus (III) - m < 0
        for λ=l:C.L
            for μ = -(λ-l):-1
                p = rlm(λ-l,μ,x,y,z)
                sum += ω₋(C[λ,μ],λ,μ,l,0,1)*p((x,y,z)=>(vx,vy,vz))
            end
        end
        cShifted[l,0] = sqrt(4*pi/(2*l+1))*sum
        sum = 0
    end

    # cShifted[l,m] for l ≠ 0, m > 0
    for l=1:C.L
        for m=1:l
            # 1: Summand 1 aus (I) - m > 0
            for λ=l:C.L
                for μ = m:m-(l-λ)
                    p = rlm(λ-l,μ-m,x,y,z)
                    sum += ω₊(C[λ,μ],λ,μ,l,m,1)*p((x,y,z)=>(vx,vy,vz))
                end
            end
            # 2: Summand 2 aus (I) - m > 0
            for λ=l+1:C.L
                for μ = max(1,m-(λ-l)):m-1
                    p = rlm(λ-l,abs(μ-m),x,y,z)
                    sum += ω₊(C[λ,μ],λ,μ,l,m,2)*p((x,y,z)=>(vx,vy,vz))
                end
            end
            # 3: Summand 3 aus (I) - m > 0
            for λ=l+1:C.L
                for μ = 1:-m-(l-λ)
                    p = rlm(λ-l,μ+m,x,y,z)
                    sum += ω₊(C[λ,μ],λ,μ,l,-m,3)*p((x,y,z)=>(vx,vy,vz))
                end
            end
            # 4: Summand aus (II) - m = 0
            for λ=m+l:C.L
                p = rlm(λ-l,m,x,y,z)
                sum += (ω₀(C[λ,0],λ,l,m,2)+ω₀(C[λ,0],λ,l,-m,3))*p((x,y,z)=>(vx,vy,vz))
            end
            # 5: Summand 1 aus (III) - m < 0
            for λ=l:C.L
                for μ = -m-(λ-l):-m-1 #min(-1,-m)
                    p = rlm(λ-l,μ+m,x,y,z)
                    sum += ω₋(C[λ,μ],λ,μ,l,-m,1)*p((x,y,z)=>(vx,vy,vz))
                end
            end
            # 6: Summand 2 aus (III) - m < 0
            for λ=l+1:C.L
                for μ = -m+1:min(-1,-m-(l-λ))
                    p = rlm(λ-l,-(μ+m),x,y,z)
                    sum -= ω₋(C[λ,μ],λ,μ,l,-m,2)*p((x,y,z)=>(vx,vy,vz))
                end
            end
            # 7: Summand 3 aus (III) - m < 0
            for λ=l+1:C.L
                for μ = m-(λ-l):-1
                    p = rlm(λ-l,μ-m,x,y,z)
                    sum += ω₋(C[λ,μ],λ,μ,l,m,3)*p((x,y,z)=>(vx,vy,vz))
                end
            end
            cShifted[l,m] = sqrt(4*pi/(2*l+1))*sum
            sum = 0
        end
    end

    # cShifted[l,m] for l ≠ 0, m < 0
    for l=1:C.L
        for m=-l:-1
            # 1: Summand 1 aus (I) - m > 0
            for λ=l:C.L
                for μ = -m+1:-m-(l-λ)
                    p = rlm(λ-l,-(μ+m),x,y,z)
                    sum -= ω₊(C[λ,μ],λ,μ,l,-m,1)*p((x,y,z)=>(vx,vy,vz))
                end
            end
            # 2: Summand 2 aus (I) - m > 0
            for λ=l+1:C.L
                for μ = max(1,-m-(λ-l)):-m-1
                    p = rlm(λ-l,-abs(μ+m),x,y,z)
                    sum += ω₊(C[λ,μ],λ,μ,l,-m,2)*p((x,y,z)=>(vx,vy,vz))
                end
            end
            # 3: Summand 3 aus (I) - m > 0
            for λ=l+1:C.L
                for μ = 1:m-(l-λ)
                    p = rlm(λ-l,-(μ-m),x,y,z)
                    sum += ω₊(C[λ,μ],λ,μ,l,m,3)*p((x,y,z)=>(vx,vy,vz))
                end
            end
            # 4: Summand aus (II) - m = 0
            for λ=l-m:C.L
                p = rlm(λ-l,m,x,y,z)
                sum += (ω₀(C[λ,0],λ,l,-m,2)+ω₀(C[λ,0],λ,l,m,3))*p((x,y,z)=>(vx,vy,vz))
            end
            # 5: Summand 1 aus (III) - m < 0
            for λ=l:C.L
                for μ = m-(λ-l):m#-1 #min(-1,m)
                    p = rlm(λ-l,abs(μ-m),x,y,z)
                    sum += ω₋(C[λ,μ],λ,μ,l,m,1)*p((x,y,z)=>(vx,vy,vz))
                end
            end
            # 6: Summand 2 aus (III) - m < 0
            for λ=l+1:C.L
                for μ = m+1:min(-1,m-(l-λ))
                    p = rlm(λ-l,μ-m,x,y,z)
                    sum += ω₋(C[λ,μ],λ,μ,l,m,2)*p((x,y,z)=>(vx,vy,vz))
                end
            end
            # 7: Summand 3 aus (III) - m < 0
            for λ=l+1:C.L
                for μ = -m-(λ-l):-1
                    p = rlm(λ-l,-(μ+m),x,y,z)
                    sum -= ω₋(C[λ,μ],λ,μ,l,-m,3)*p((x,y,z)=>(vx,vy,vz))
                end
            end
            cShifted[l,m] = sqrt(4*pi/(2*l+1))*sum
            sum = 0
        end
    end

    # turn spherical coefficients into solid coefficients
    if cShifted.solid
        cShifted.solid = false
        solid!(cShifted)
        solid!(C)
    end

    return cShifted
end

function ω₀(c,λ,l,m,num)
    ω = c*factorial(λ)*sqrt((2*λ+1)/(4*pi))/sqrt(4*factorial(l+m)*factorial(l-m)*factorial(λ-l+m)*factorial(λ-l-m))
    if num == 1 # m == 0
        ω *= 2
    else # num == 2 || num == 3
        ω *= (-1)^m
    end
    return ω
end

function ω₊(c,λ,μ,l,m,num)
    ω = c*sqrt(2*factorial(λ+μ)*factorial(λ-μ))*sqrt((2*λ+1)/(4*pi))/sqrt(4*factorial(l+m)*factorial(l-m)*factorial(λ-l+μ-m)*factorial(λ-l-(μ-m)))
    if num == 1
        if m == μ || m == 0
            ω *= sqrt(2)
        else
            ω *= 1
        end
    elseif num == 2
        ω *= (-1)^(μ-m)
    else # num == 3
        ω *= (-1)^m
    end
    return ω
end

function ω₋(c,λ,μ,l,m,num)
    ω = c*(-1)^μ*sqrt(2*factorial(λ+μ)*factorial(λ-μ))*sqrt((2*λ+1)/(4*pi))/sqrt(4*factorial(l+m)*factorial(l-m)*factorial(λ-l+μ-m)*factorial(λ-l-(μ-m)))
    if num == 1
        if m == μ || m == 0
            ω *= (-1)^μ*sqrt(2)
        else
            ω *= (-1)^μ
        end
    elseif num == 2
        ω *= (-1)^m
    else # num == 3
        ω *= (-1)^(μ-m)
    end
    return ω
end


function translateRlm(l::Int64, m::Int64,vx,vy,vz)

    @polyvar x y z

    Rlmt = 0;

    # two cases m≥0 and m<0, each with three summands:
    if m >= 0
    for λ=0:l
        for μ=max(0,λ-l+m):min(λ,m)
            if m == 0
                p = rlm(l-λ,0,x,y,z);
                mult = 1/(factorial(λ)*factorial(l-λ));
                sum1 = rlm(λ,0,x,y,z)*p((x,y,z)=>(vx,vy,vz));
                Rlmt += mult*sum1;
            elseif μ == 0 # => m-μ = m > 0
                p = rlm(l-λ,m,x,y,z);
                mult = 1/(factorial(λ)*sqrt(2*factorial(l-λ+m)*factorial(l-λ-m)));
                sum1 = rlm(λ,0,x,y,z)*p((x,y,z)=>(vx,vy,vz));
                Rlmt += mult*sum1;
            elseif μ == m # => m-μ = 0
                p = rlm(l-λ,0,x,y,z);
                mult = 1/(factorial(l-λ)*sqrt(2*factorial(λ+m)factorial(λ-m)));
                sum1 = rlm(λ,m,x,y,z)*p((x,y,z)=>(vx,vy,vz));
                Rlmt += mult*sum1;
            else # m-μ > 0, μ > 0
                p = rlm(l-λ,m-μ,x,y,z);
                q = rlm(l-λ,-(m-μ),x,y,z);
                mult = 1/sqrt(4*factorial(λ+μ)*factorial(λ-μ)*factorial(l-λ+m-μ)*factorial(l-λ-(m-μ)));
                sum1 = rlm(λ,μ,x,y,z)*p((x,y,z)=>(vx,vy,vz));
                sum2 = rlm(λ,-μ,x,y,z)*q((x,y,z)=>(vx,vy,vz));
                Rlmt += mult*(sum1 - sum2);
            end
        end
    end

    for λ=m+1:l-1
        for μ=m+1:min(λ,-λ+l+m) # μ > 0, m-μ < 0
            p = rlm(l-λ,abs(m-μ),x,y,z);
            q = rlm(l-λ,-abs(m-μ),x,y,z);
            mult = (-1)^(m-μ)/sqrt(4*factorial(λ+μ)*factorial(λ-μ)*factorial(l-λ+(m-μ))*factorial(l-λ-(m-μ)));
            sum1 = rlm(λ,μ,x,y,z)*p((x,y,z)=>(vx,vy,vz));
            sum2 = rlm(λ,-μ,x,y,z)*q((x,y,z)=>(vx,vy,vz));
            Rlmt += mult*(sum1 + sum2);
        end
    end

    for λ=1:l-m-1
        for μ=max(-λ,λ-l+m):-1 # μ < 0, m-μ > 0
            p = rlm(l-λ,m-μ,x,y,z);
            q = rlm(l-λ,-(m-μ),x,y,z);
            mult = (-1)^(μ)/sqrt(4*factorial(λ+μ)*factorial(λ-μ)*factorial(l-λ+m-μ)*factorial(l-λ-(m-μ)));
            sum1 = rlm(λ,abs(μ),x,y,z)*p((x,y,z)=>(vx,vy,vz));
            sum2 = rlm(λ,-abs(μ),x,y,z)*q((x,y,z)=>(vx,vy,vz));
            Rlmt += mult*(sum1 + sum2);
        end
    end

    Rlmt *= sqrt(factorial(l-m)*factorial(l+m));
    Rlmt *= (m == 0) ? 1 : sqrt(2);

    else # m < 0
        for λ=0:l
            for μ=max(m,-λ):min(0,-λ+l+m)
                if μ == 0 # => m-μ = m < 0
                    p = rlm(l-λ,-abs(m),x,y,z);
                    mult = (-1)^m/(factorial(λ)*sqrt(2*factorial(l-λ+m)*factorial(l-λ-m)));
                    sum1 = rlm(λ,0,x,y,z)*p((x,y,z)=>(vx,vy,vz));
                    Rlmt += (-1)*mult*sum1;
                elseif μ == m # => m-μ = 0
                    p = rlm(l-λ,0,x,y,z);
                    mult = (-1)^m/(factorial(l-λ)*sqrt(2*factorial(λ+m)factorial(λ-m)));
                    sum1 = rlm(λ,-abs(m),x,y,z)*p((x,y,z)=>(vx,vy,vz));
                    Rlmt += (-1)*mult*sum1;
                else # μ < 0, m-μ < 0
                    p = rlm(l-λ,abs(m-μ),x,y,z);
                    q = rlm(l-λ,-abs(m-μ),x,y,z);
                    mult = (-1)^m/sqrt(4*factorial(λ+μ)*factorial(λ-μ)*factorial(l-λ+(m-μ))*factorial(l-λ-(m-μ)));
                    sum1 = rlm(λ,-abs(μ),x,y,z)*p((x,y,z)=>(vx,vy,vz));
                    sum2 = rlm(λ,abs(μ),x,y,z)*q((x,y,z)=>(vx,vy,vz));
                    Rlmt += (-1)*mult*(sum1 + sum2);
                end
            end
        end

        for λ=-m+1:l-1
            for μ=max(-λ,λ-l+m):m-1 # μ < 0, m-μ > 0
                p = rlm(l-λ,-(m-μ),x,y,z);
                q = rlm(l-λ,(m-μ),x,y,z);
                mult = (-1)^(μ)/sqrt(4*factorial(λ+μ)*factorial(λ-μ)*factorial(l-λ+m-μ)*factorial(l-λ-(m-μ)));
                sum1 = rlm(λ,abs(μ),x,y,z)*p((x,y,z)=>(vx,vy,vz));
                sum2 = rlm(λ,-abs(μ),x,y,z)*q((x,y,z)=>(vx,vy,vz));
                Rlmt += mult*(sum1 - sum2);
            end
        end

        for λ=1:l+m-1
            for μ=1:min(λ,-λ+l+m) # μ > 0, m-μ < 0
                p = rlm(l-λ,abs(m-μ),x,y,z);
                q = rlm(l-λ,-abs(m-μ),x,y,z);
                mult = (-1)^(m-μ)/sqrt(4*factorial(λ+μ)*factorial(λ-μ)*factorial(l-λ+(m-μ))*factorial(l-λ-(m-μ)));
                sum1 = rlm(λ,-μ,x,y,z)*p((x,y,z)=>(vx,vy,vz));
                sum2 = rlm(λ,μ,x,y,z)*q((x,y,z)=>(vx,vy,vz));
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

    cShifted = deepcopy(C)

    vx = v[1]
    vy = v[2]
    vz = v[3]

    # turn solid coefficients into spherical coefficients
    if C.solid
        spherical!(C)
    end

    # cShifted[0,0]:
    sum = 0
    # Summand aus (I) - m > 0
    for λ=1:C.L
        for μ = 1:λ
            p = rlm(λ,μ,x,y,z)
            sum += abs(ω₊(C[λ,μ],λ,μ,0,0,1)*p((x,y,z)=>(vx,vy,vz)))
        end
    end
    # Summand aus (II) - m = 0
    for λ=0:C.L
        p = rlm(λ,0,x,y,z)
        sum += abs(ω₀(C[λ,0],λ,0,0,1)*p((x,y,z)=>(vx,vy,vz)))
    end
    # Summand aus (III) - m < 0
    for λ=1:C.L
        for μ = -λ:-1
            p = rlm(λ,μ,x,y,z)
            sum += abs(ω₋(C[λ,μ],λ,μ,0,0,1)*p((x,y,z)=>(vx,vy,vz)))
        end
    end
    cShifted[0,0] = sqrt(4*pi)*sum
    sum = 0

    # cShifted[l,0] (l > 0):
    for l=1:C.L
        # Summand aus (I) - m > 0
        for λ=l:C.L
            for μ = 1:λ-l
                p = rlm(λ-l,μ,x,y,z)
                sum += abs(ω₊(C[λ,μ],λ,μ,l,0,1)*p((x,y,z)=>(vx,vy,vz)))
            end
        end
        # Summand aus (II) - m = 0
        for λ=l:C.L
            p = rlm(λ-l,0,x,y,z)
            sum += abs(ω₀(C[λ,0],λ,l,0,1)*p((x,y,z)=>(vx,vy,vz)))
        end
        # Summand aus (III) - m < 0
        for λ=l:C.L
            for μ = -(λ-l):-1
                p = rlm(λ-l,μ,x,y,z)
                sum += abs(ω₋(C[λ,μ],λ,μ,l,0,1)*p((x,y,z)=>(vx,vy,vz)))
            end
        end
        cShifted[l,0] = sqrt(4*pi/(2*l+1))*sum
        sum = 0
    end

    # cShifted[l,m] for l ≠ 0, m > 0
    for l=1:C.L
        for m=1:l
            # 1: Summand 1 aus (I) - m > 0
            for λ=l:C.L
                for μ = m:m-(l-λ)
                    p = rlm(λ-l,μ-m,x,y,z)
                    sum += abs(ω₊(C[λ,μ],λ,μ,l,m,1)*p((x,y,z)=>(vx,vy,vz)))
                end
            end
            # 2: Summand 2 aus (I) - m > 0
            for λ=l+1:C.L
                for μ = max(1,m-(λ-l)):m-1
                    p = rlm(λ-l,abs(μ-m),x,y,z)
                    sum += abs(ω₊(C[λ,μ],λ,μ,l,m,2)*p((x,y,z)=>(vx,vy,vz)))
                end
            end
            # 3: Summand 3 aus (I) - m > 0
            for λ=l+1:C.L
                for μ = 1:-m-(l-λ)
                    p = rlm(λ-l,μ+m,x,y,z)
                    sum += abs(ω₊(C[λ,μ],λ,μ,l,-m,3)*p((x,y,z)=>(vx,vy,vz)))
                end
            end
            # 4: Summand aus (II) - m = 0
            for λ=m+l:C.L
                p = rlm(λ-l,m,x,y,z)
                sum += abs((ω₀(C[λ,0],λ,l,m,2)+ω₀(C[λ,0],λ,l,-m,3))*p((x,y,z)=>(vx,vy,vz)))
            end
            # 5: Summand 1 aus (III) - m < 0
            for λ=l:C.L
                for μ = -m-(λ-l):-m-1 #min(-1,-m)
                    p = rlm(λ-l,μ+m,x,y,z)
                    sum += abs(ω₋(C[λ,μ],λ,μ,l,-m,1)*p((x,y,z)=>(vx,vy,vz)))
                end
            end
            # 6: Summand 2 aus (III) - m < 0
            for λ=l+1:C.L
                for μ = -m+1:min(-1,-m-(l-λ))
                    p = rlm(λ-l,-(μ+m),x,y,z)
                    sum += abs(ω₋(C[λ,μ],λ,μ,l,-m,2)*p((x,y,z)=>(vx,vy,vz)))
                end
            end
            # 7: Summand 3 aus (III) - m < 0
            for λ=l+1:C.L
                for μ = m-(λ-l):-1
                    p = rlm(λ-l,μ-m,x,y,z)
                    sum += abs(ω₋(C[λ,μ],λ,μ,l,m,3)*p((x,y,z)=>(vx,vy,vz)))
                end
            end
            cShifted[l,m] = sqrt(4*pi/(2*l+1))*sum
            sum = 0
        end
    end

    # cShifted[l,m] for l ≠ 0, m < 0
    for l=1:C.L
        for m=-l:-1
            # 1: Summand 1 aus (I) - m > 0
            for λ=l:C.L
                for μ = -m+1:-m-(l-λ)
                    p = rlm(λ-l,-(μ+m),x,y,z)
                    sum += abs(ω₊(C[λ,μ],λ,μ,l,-m,1)*p((x,y,z)=>(vx,vy,vz)))
                end
            end
            # 2: Summand 2 aus (I) - m > 0
            for λ=l+1:C.L
                for μ = max(1,-m-(λ-l)):-m-1
                    p = rlm(λ-l,-abs(μ+m),x,y,z)
                    sum += abs(ω₊(C[λ,μ],λ,μ,l,-m,2)*p((x,y,z)=>(vx,vy,vz)))
                end
            end
            # 3: Summand 3 aus (I) - m > 0
            for λ=l+1:C.L
                for μ = 1:m-(l-λ)
                    p = rlm(λ-l,-(μ-m),x,y,z)
                    sum += abs(ω₊(C[λ,μ],λ,μ,l,m,3)*p((x,y,z)=>(vx,vy,vz)))
                end
            end
            # 4: Summand aus (II) - m = 0
            for λ=l-m:C.L
                p = rlm(λ-l,m,x,y,z)
                sum += abs((ω₀(C[λ,0],λ,l,-m,2)+ω₀(C[λ,0],λ,l,m,3))*p((x,y,z)=>(vx,vy,vz)))
            end
            # 5: Summand 1 aus (III) - m < 0
            for λ=l:C.L
                for μ = m-(λ-l):m#-1 #min(-1,m)
                    p = rlm(λ-l,abs(μ-m),x,y,z)
                    sum += abs(ω₋(C[λ,μ],λ,μ,l,m,1)*p((x,y,z)=>(vx,vy,vz)))
                end
            end
            # 6: Summand 2 aus (III) - m < 0
            for λ=l+1:C.L
                for μ = m+1:min(-1,m-(l-λ))
                    p = rlm(λ-l,μ-m,x,y,z)
                    sum += abs(ω₋(C[λ,μ],λ,μ,l,m,2)*p((x,y,z)=>(vx,vy,vz)))
                end
            end
            # 7: Summand 3 aus (III) - m < 0
            for λ=l+1:C.L
                for μ = -m-(λ-l):-1
                    p = rlm(λ-l,-(μ+m),x,y,z)
                    sum += abs(ω₋(C[λ,μ],λ,μ,l,-m,3)*p((x,y,z)=>(vx,vy,vz)))
                end
            end
            cShifted[l,m] = sqrt(4*pi/(2*l+1))*sum
            sum = 0
        end
    end

    # turn spherical coefficients into solid coefficients
    if cShifted.solid
        cShifted.solid = false
        solid!(cShifted)
        solid!(C)
    end

    return cShifted
end
