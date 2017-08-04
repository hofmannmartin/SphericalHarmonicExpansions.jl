@testset "spherical harmonics expansion" begin
  ɛ = eps(Float64)
  @polyvar x y z

  #Test by using unit vectors
  #l = 10 -> size(C) = l²+2l+1 = 121
  for l in 0:10
    for m in -l:l
      C = zeros(121);
      C[l*(l+1)+m+1] = 1;
      @test isapprox(sphericalHarmonicsExpansion(C,x,y,z),ylm(l,m,x,y,z),atol=ɛ)
    end
  end

  #Testing the addition
  C = [1;0.5;10.78;0.87456]
  @test isapprox(sphericalHarmonicsExpansion(C,x,y,z),1*ylm(0,0,x,y,z)+0.5*ylm(1,-1,x,y,z)+10.78*ylm(1,0,x,y,z)+0.87456*ylm(1,1,x,y,z),atol=ɛ)
end
