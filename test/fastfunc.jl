@testset "macro for generating fast functions" begin
	ε = 100*eps(Float64)
	f(x1,x2,x3) = 3*x3^2 + 2*x2*x1
	@polyvar x y z
	polynomial = 3*z^2 + 2*y*x
	g = @fastfunc polynomial

	t1,t2,t3 = randn(3)
	@test isapprox(f(t1,t2,t3),g(t1,t2,t3),atol=ɛ)
end
