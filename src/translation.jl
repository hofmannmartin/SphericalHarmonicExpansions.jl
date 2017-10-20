# Translation of the coefficients c by a vector v
function translation(C::SphericalHarmonicCoefficients,v,x::Variable, y::Variable, z::Variable)

    c = SphericalHarmonicCoefficients(zeros((C.L+1)^2))

    vx = v[1]
    vy = v[2]
    vz = v[3]

    # c₀₀:
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
    c[0,0] = sqrt(4*pi)*sum
    sum = 0

    # cₗ₀ (l > 0):
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
        c[l,0] = sqrt(4*pi/(2*l+1))*sum
        sum = 0
    end

    # cₗₘ für l ≠ 0, m > 0
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
            c[l,m] = sqrt(4*pi/(2*l+1))*sum
            sum = 0
        end
    end

    # cₗₘ für l ≠ 0, m < 0
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
            c[l,m] = sqrt(4*pi/(2*l+1))*sum
            sum = 0
        end
    end

    return c
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


function translateRlm(l::Int64, m::Int64,vx,vy,vz, x::Variable, y::Variable, z::Variable)

    Rlmt = 0;

    # 3 Summanden fuer m≥0 und m<0 getrennt:
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
                    p = rlm(l-λ,m,x,y,z);
                    mult = (-1)^m/(factorial(λ)*sqrt(2*factorial(l-λ+m)*factorial(l-λ-m)));
                    sum1 = rlm(λ,0,x,y,z)*p((x,y,z)=>(vx,vy,vz));
                    Rlmt += mult*sum1;
                elseif μ == m # => m-μ = 0
                    p = rlm(l-λ,0,x,y,z);
                    mult = (-1)^m/(factorial(l-λ)*sqrt(2*factorial(λ+m)factorial(λ-m)));
                    sum1 = rlm(λ,m,x,y,z)*p((x,y,z)=>(vx,vy,vz));
                    Rlmt += mult*sum1;
                else # μ < 0, m-μ < 0
                    p = rlm(l-λ,(m-μ),x,y,z);
                    q = rlm(l-λ,-(m-μ),x,y,z);
                    mult = (-1)^m/sqrt(4*factorial(λ+μ)*factorial(λ-μ)*factorial(l-λ+(m-μ))*factorial(l-λ-(m-μ)));
                    sum1 = rlm(λ,abs(μ),x,y,z)*p((x,y,z)=>(vx,vy,vz));
                    sum2 = rlm(λ,-abs(μ),x,y,z)*q((x,y,z)=>(vx,vy,vz));
                    Rlmt += mult*(sum1 + sum2);
                end
            end
        end

        for λ=-m+1:l-1
            for μ=max(-λ,λ-l+m):m-1 # μ < 0, m-μ > 0
                p = rlm(l-λ,m-μ,x,y,z);
                q = rlm(l-λ,-(m-μ),x,y,z);
                mult = (-1)^(μ)/sqrt(4*factorial(λ+μ)*factorial(λ-μ)*factorial(l-λ+m-μ)*factorial(l-λ-(m-μ)));
                sum1 = rlm(λ,μ,x,y,z)*p((x,y,z)=>(vx,vy,vz));
                sum2 = rlm(λ,-μ,x,y,z)*q((x,y,z)=>(vx,vy,vz));
                Rlmt += mult*(sum1 - sum2);
            end
        end

        for λ=1:l+m-1
            for μ=1:min(λ,-λ+l+m) # μ > 0, m-μ < 0
                p = rlm(l-λ,(m-μ),x,y,z);
                q = rlm(l-λ,-(m-μ),x,y,z);
                mult = (-1)^(m-μ)/sqrt(4*factorial(λ+μ)*factorial(λ-μ)*factorial(l-λ+(m-μ))*factorial(l-λ-(m-μ)));
                sum1 = rlm(λ,μ,x,y,z)*p((x,y,z)=>(vx,vy,vz));
                sum2 = rlm(λ,-μ,x,y,z)*q((x,y,z)=>(vx,vy,vz));
                Rlmt += mult*(sum1 - sum2);
            end
        end

        Rlmt *= (-1)^m*sqrt(2*factorial(l-abs(m))*factorial(l+abs(m)));

    # else # m = 0
    #     for λ=0:l
    #         # μ = 0
    #         p = rlm(l-λ,0,x,y,z);
    #         mult = 1/(factorial(λ)*factorial(l-λ));
    #         sum1 = rlm(λ,0,x,y,z)*p((x,y,z)=>(vx,vy,vz));
    #         Rlmt += mult*sum1;
    #     end
    #
    #     for λ=1:l-1
    #         for μ=1:min(λ,-λ+l) # μ > 0, m-μ = -μ < 0
    #             p = rlm(l-λ,μ,x,y,z);
    #             q = rlm(l-λ,-μ,x,y,z);
    #             mult = (-1)^(μ)/sqrt(4*factorial(λ+μ)*factorial(λ-μ)*factorial(l-λ+μ)*factorial(l-λ-μ));
    #             sum1 = rlm(λ,μ,x,y,z)*p((x,y,z)=>(vx,vy,vz));
    #             sum2 = rlm(λ,-μ,x,y,z)*q((x,y,z)=>(vx,vy,vz));
    #             Rlmt += mult*(sum1 + sum2);
    #         end
    #
    #         for μ=max(-λ,λ-l):-1 # μ < 0, m-μ = -μ > 0
    #             p = rlm(l-λ,-μ,x,y,z);
    #             q = rlm(l-λ,μ,x,y,z);
    #             mult = (-1)^(μ)/sqrt(4*factorial(λ-μ)*factorial(λ+μ)*factorial(l-λ-μ)*factorial(l-λ+μ));
    #             sum1 = rlm(λ,-μ,x,y,z)*p((x,y,z)=>(vx,vy,vz));
    #             sum2 = rlm(λ,μ,x,y,z)*q((x,y,z)=>(vx,vy,vz));
    #             Rlmt += mult*(sum1 + sum2);
    #         end
    #     end
    #
    #     Rlmt *= factorial(l);
    end
    return Rlmt
end
