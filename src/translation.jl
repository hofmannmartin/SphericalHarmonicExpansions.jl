# adding r^l to a spherical harmonic and canceeling out with the normalization of the coordinates
function addRl(l::Int64, ylm::Polynomial, r::Variable)
    # Zerlegung des Polynoms in Terme:
    Summanden = terms(ylm);

    rylm = 0;

    for i=1:length(Summanden)
        deg = degree(monomial(Summanden[i])); # Gibt den gesamten Grad des Monoms an
        degR = l-deg; # durch das Kürzen ergibt sich ein Grad von l-deg fuer r
        Summanden[i] = r^degR*Summanden[i]; # Terme von ylm mit Faktor r
        rylm += Summanden[i]; # Addition der Terme zu einem Polynom
    end

    return rylm;
end

# solid harmonics
function rlm(l::Int64, m::Int64, r::Variable, x::Variable, y::Variable, z::Variable)
    rlm = addRl(l,ylm(l,m,x,y,z)+0*r,r);
    rlm = sqrt(4*pi/(2*l+1))*rlm;
    #rlm = sqrt(4*pi/(2*l+1))*r^l*ylm(l,m,x,y,z); # (x,y,z) sind normiert, r Normierung
    return rlm;
end

# Translation of the coefficients c by a vector v
function translation(C::SphericalHarmonicCoefficients,vx,vy,vz)
    @polyvar r x y z;
    Ct = SphericalHarmonicCoefficients(zeros((C.L+1)^2));
    # Normierung von v:
        vr = sqrt(vx^2+vy^2+vz^2);
        vxn = vx/vr;
        vyn = vy/vr;
        vzn = vz/vr;
    for l=0:C.L
        for m=-l:l
            for λ=l:C.L
                sum = 0;
                for μ=-λ:λ
                    if binomial(λ+μ,l+m) != 0 && binomial(λ-μ,l-m) != 0 # else: rlm(λ-l,μ-m,x,y,z) undefined
                        p = rlm(λ-l,μ-m,r,x,y,z);
                        sum += C[λ,μ]*p((r,x,y,z)=>(vr,vxn,vyn,vzn))*sqrt(binomial(λ+μ,l+m)*binomial(λ-μ,l-m));
                    end
                end
                Ct[l,m] += sqrt((2*λ+1)/(2*l+1))*sum;
            end
        end
    end
    return Ct;
end

function translateRlm(l::Int64, m::Int64,vx,vy,vz, r::Variable, x::Variable, y::Variable, z::Variable)
#    @polyvar x y z;
    Rlmt = 0;
    # Normierung von v:
        vr = sqrt(vx^2+vy^2+vz^2);
        vxn = vx/vr;
        vyn = vy/vr;
        vzn = vz/vr;

    for λ=0:l
        for μ=-λ:λ
                if abs(m-μ) <= l-λ
                    p = rlm(l-λ,m-μ,r,x,y,z);
                else
                    p = 0*x;
                end
                Rlmt += rlm(λ,μ,r,x,y,z)*p((r,x,y,z)=>(vr,vxn,vyn,vzn))*sqrt(binomial(l+m,λ+μ)*binomial(l-m,λ-μ));
        end
    end
    return Rlmt;
end
