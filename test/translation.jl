@testset "multiply with rˡ" begin
ε = 32*eps(Float64);
@polyvar r x y z;

@test isapprox(addRl(3,1+0*r,r),r^3+0,atol=ε);
@test isapprox(addRl(7,0*r+0,r),0,atol=ε);
@test isapprox(addRl(5,1+2*x+3*x*y+4*z^3+5*x*y^2*z+6*x^3*z^2+0*r,r),r^5+2*x*r^4+3*x*y*r^3+4*z^3*r^2+5*x*y^2*z*r+6*x^3*z^2,atol=ε);
end

@testset "solid harmonic" begin
ε = 32*eps(Float64);
@polyvar r x y z;

# test for l = 0,1,2,3

# l = 0
@test isapprox(rlm(0,0,r,x,y,z),1+0*x,atol=ε);

# l = 1
@test isapprox(rlm(1,-1,r,x,y,z),y+0,atol=ε);
@test isapprox(rlm(1,0,r,x,y,z),z+0,atol=ε);
@test isapprox(rlm(1,1,r,x,y,z),x+0,atol=ε);

# l = 2
@test isapprox(rlm(2,-2,r,x,y,z),sqrt(3)*x*y+0,atol=ε);
@test isapprox(rlm(2,-1,r,x,y,z),sqrt(3)*y*z+0,atol=ε);
@test isapprox(rlm(2,0,r,x,y,z),1/2*(3*z^2-r^2),atol=ε);
@test isapprox(rlm(2,1,r,x,y,z),sqrt(3)*x*z+0,atol=ε);
@test isapprox(rlm(2,2,r,x,y,z),sqrt(3/4)*(x^2-y^2),atol=ε);

# l = 3
@test isapprox(rlm(3,-3,r,x,y,z),sqrt(5/8)*(3*x^2*y-y^3),atol=ε);
@test isapprox(rlm(3,-2,r,x,y,z),sqrt(15)*x*y*z+0,atol=ε);
@test isapprox(rlm(3,-1,r,x,y,z),sqrt(3/8)*y*(5*z^2-r^2),atol=ε);
@test isapprox(rlm(3,0,r,x,y,z),1/2*(5*z^3-3*z*r^2),atol=ε);
@test isapprox(rlm(3,1,r,x,y,z),sqrt(3/8)*x*(5*z^2-r^2),atol=ε);
@test isapprox(rlm(3,2,r,x,y,z),sqrt(15/4)*(x^2-y^2)*z,atol=ε);
@test isapprox(rlm(3,3,r,x,y,z),sqrt(5/8)*(x^3-3*x*y^2),atol=ε);

end
