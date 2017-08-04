@testset "spherical harmonics expansion" begin
  ɛ = eps(Float64)
  @polyvar x y z

  @test_throws DomainError SphericalHarmonicCoefficients(-2)
  @test_throws DomainError SphericalHarmonicCoefficients([1;2;3;4;5])

  c1 = SphericalHarmonicCoefficients(2)
  c2 = SphericalHarmonicCoefficients([1;2;3;4])

  @test isapprox(length(c1.c),9,atol=ɛ)
  @test isapprox(c2.L, 1)

  @test isapprox(c2[1],1,atol=ɛ)
  @test isapprox(c2[2],2,atol=ɛ)
  @test isapprox(c2[3],3,atol=ɛ)
  @test isapprox(c2[4],4,atol=ɛ)

  @test isapprox(c2[0,0],1,atol=ɛ)
  @test isapprox(c2[1,-1],2,atol=ɛ)
  @test isapprox(c2[1,0],3,atol=ɛ)
  @test isapprox(c2[1,1],4,atol=ɛ)

  c2[1] = 10
  c2[2] = 20
  c2[3] = 30
  c2[4] = 40
  @test isapprox(c2[1],10,atol=ɛ)
  @test isapprox(c2[2],20,atol=ɛ)
  @test isapprox(c2[3],30,atol=ɛ)
  @test isapprox(c2[4],40,atol=ɛ)

  c2[0,0] = 100
  c2[1,-1] = 200
  c2[1,0] = 300
  c2[1,1] = 400
  @test isapprox(c2[1],100,atol=ɛ)
  @test isapprox(c2[2],200,atol=ɛ)
  @test isapprox(c2[3],300,atol=ɛ)
  @test isapprox(c2[4],400,atol=ɛ)

  #Test by using unit vectors
  #l = 10 -> size(C) = l²+2l+1 = 121

  for l in 0:2
    for m in -l:l
      C = SphericalHarmonicCoefficients(2);
      C[l,m] = 1;
      @test isapprox(sphericalHarmonicsExpansion(C,x,y,z),ylm(l,m,x,y,z),atol=ɛ)
    end
  end

  #Testing the addition
  C = SphericalHarmonicCoefficients([1;0.5;10.78;0.87456])
  @test isapprox(sphericalHarmonicsExpansion(C,x,y,z),1*ylm(0,0,x,y,z)+0.5*ylm(1,-1,x,y,z)+10.78*ylm(1,0,x,y,z)+0.87456*ylm(1,1,x,y,z),atol=ɛ)
end
