@testset "Translation of Rlm" begin

ε = 1000*eps(Float64)
@polyvar x y z;

Wert = zeros(25,2)
v = randn(3) # Verschiebung
pkt = randn(3) # Testpunkt

vp = pkt+v;
              
for l=0:4
    for m=-l:l
        q = rlm(l,m,x,y,z);
        p = translateRlm(l,m,v[1],v[2],v[3],x,y,z);

        polyRlm = @fastfunc q;
        polyTranslateRlm = @fastfunc p;

        Wert[l*(l+1)+m+1,1] = polyRlm(vp[1],vp[2],vp[3]);
        Wert[l*(l+1)+m+1,2] = polyTranslateRlm(pkt[1],pkt[2],pkt[3]);
    end
end

@test isapprox(Wert[:,1],Wert[:,2],atol = ε)
end

@testset "Translation of Clm" begin

@polyvar x y z

v = randn(3) # Verschiebung

L = 6
C = SphericalHarmonicCoefficients(randn((L+1)^2))

# Translation der Koeffizienten mittels v, danach zurück mit -v
coeffsTrans = translation(C,v,x,y,z)
coeffsBack = translation(coeffsTrans,-v,x,y,z)

@test isapprox(C,coeffsBack)

# Konkretes Beispiel für die Funktion xyz
C = SphericalHarmonicCoefficients(zeros(25)); 
C[3,-2] = 1/sqrt(15)*sqrt((4*pi)/(2*3+1));

v = [1,-1,2]
coeffsTrans = translation(C,v,x,y,z)

Ctrans = zeros(25)
Ctrans[1:11] = sqrt(4pi/3).*[-2*sqrt(3), 2, -1, -2, 2/sqrt(5), 1/sqrt(5), 0, -1/sqrt(5), 0, 0, 1/sqrt(35)]
Ctrans = SphericalHarmonicCoefficients(Ctrans)

@test isapprox(Ctrans,coeffsTrans)

end
