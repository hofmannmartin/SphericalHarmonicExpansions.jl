@testset "Generating fast functions" begin
	ε = 100*eps(Float64)
	r = randn(3)
	t1,t2,t3 = r

	f(x1,x2,x3) = 3*x3^2 + 2*x2*x1
	@polyvar x y z
	polynomial = 3*z^2 + 2*y*x
	g = @fastfunc polynomial
	@test isapprox(f(t1,t2,t3),g(t1,t2,t3),atol=ɛ)
	# Test inside function scope
	function useInsideFunctionScope1(polynomial,t1,t2,t3)
		g = @fastfunc polynomial
		@test isapprox(f(t1,t2,t3),Base.invokelatest(g, t1,t2,t3),atol=ɛ)
	end
	useInsideFunctionScope1(polynomial,t1,t2,t3)

	h = fastfunc(polynomial)
    @test isapprox(f(t1,t2,t3),h(r),atol=ɛ)
	# Test inside function scope
	function useInsideFunctionScope2(polynomial,t1,t2,t3,r)
        g = fastfunc(polynomial)
		@test isapprox(f(t1,t2,t3),g(r),atol=ɛ)
	end
    useInsideFunctionScope2(polynomial,t1,t2,t3,r)

    # Test fastfunc on SphericalHarmonicCoefficients
    ca = rand(4)
    coeffs = SphericalHarmonicCoefficients(ca,1.0,true)
    fc = fastfunc(coeffs)
    @test isapprox(fc(ones(3)),sum(ca),atol=ε)
end
