##################################################
## Test of the translation of a solid expansion ##
##################################################
@testset "Translation of Rlm" begin

    ε = 1000*eps(Float64)
    @polyvar x y z;

    Wert = zeros(25,2)
    v = randn(3) # Verschiebung
    pkt = randn(3) # Testpunkt

    vp = pkt+v;
                  
    for l=0:4
        for m=-l:l
            q = SphericalHarmonics.zlm(l,m,x,y,z);
            p = SphericalHarmonics.translateRlm(l,m,v[1],v[2],v[3]);

            polyRlm = @fastfunc q;
            polyTranslateRlm = @fastfunc p;

            Wert[l*(l+1)+m+1,1] = polyRlm(vp[1],vp[2],vp[3]);
            Wert[l*(l+1)+m+1,2] = polyTranslateRlm(pkt[1],pkt[2],pkt[3]);
        end
    end

    @test isapprox(Wert[:,1],Wert[:,2],atol = ε)
end

#################################################
## Test of the translation of the coefficients ##
#################################################
@testset "Translation of Clm" begin

    @polyvar x y z

    ###################################
    # Translation of the coefficients #
    # with v and after that with -v   #
    ###################################
    # Shift vector
    v = randn(3)

    # Coefficients
    L = 6
    C = SphericalHarmonicCoefficients(randn((L+1)^2))

    # Translation
    coeffsTrans = translation(C,v)
    coeffsBack = translation(coeffsTrans,-v)

    # Test
    @test isapprox(C,coeffsBack)

    ############################################
    # Konkretes Beispiel für die Funktion xyz  #
    # solid & spherical Expansion              #
    ############################################
    # Coefficients
    Cspher = SphericalHarmonicCoefficients(zeros(25)); 
    Cspher[3,-2] = 1/sqrt(15)*sqrt((4*pi)/(2*3+1));
    Csolid = SphericalHarmonicCoefficients(zeros(25),1.0,true); 
    Csolid[3,-2] = 1/sqrt(15);

    # Translation
    v = [1.,-1.,2.]
    coeffsTransSpher = translation(Cspher,v)
    coeffsTransSolid = translation(Csolid,v)

    # Correct solution
    CtransSpher = zeros(25)
    CtransSpher[1:11] = sqrt(4pi/3).*[-2*sqrt(3), 2, -1, -2, 2/sqrt(5), 1/sqrt(5), 0, -1/sqrt(5), 0, 0, 1/sqrt(35)]
    CtransSpher = SphericalHarmonicCoefficients(CtransSpher)
    CtransSolid = zeros(25)
    CtransSolid[1:11] = [-2, 2, -1, -2, 2/sqrt(3), 1/sqrt(3), 0, -1/sqrt(3), 0, 0, 1/sqrt(15)]
    CtransSolid = SphericalHarmonicCoefficients(CtransSolid,1.0,true)

    # Test
    @test isapprox(CtransSpher,coeffsTransSpher)
    @test isapprox(CtransSolid,coeffsTransSolid)
end
