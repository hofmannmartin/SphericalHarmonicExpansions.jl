"""
    @fastfunc functionname[::String] polynomial[::Polynomial{true}]

    Defines a function with name `functionname` for `polynomial` for fast numerical evaluation.
    Note that the polynomial is transformed using `@fastmath`, which calls functions that
    may violate strict IEEE semantics.

    # Examples
    ```
    julia> using SphericalHarmonics

    julia> @polyvar x y z
    (x, y, z)

    julia> p = 15.0*x*y^2+7.5*x*z^13
    7.5xz¹³ + 15.0xy²

    julia> foo = @fastfunc p;

    julia> foo(1.0,2.0,3.0)
    1.19574825e7
    ```
    """
macro fastfunc(polynomial)
	return quote
		local polystr = string($(esc(polynomial)))
		local vars = string(tuple(variables($(esc(polynomial)))...))
        # create expression for function definition
		eval(Meta.parse(vars*" -> @fastmath "*polystr))
	end
end

function fastfunc(polynomial)
    expr = string(polynomial)
    args = string(variables(polynomial))
    return mk_function(Meta.parse(args*" -> "*expr))
end