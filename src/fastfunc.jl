"""
    @fastfunc functionname[::String] polynomial[::Polynomial{true}]

    Defines a function with name `functionname` for `polynomial` for fast numerical evaluation.
    Note that the polynomial is transformed using `@fastmath`, which calls functions that
    may violate strict IEEE semantics.

    # Examples
    ```
    julia> using MultivariatePolynomials

    julia> @polyvar x y z
    z

    julia> p = 15.0*x*y^2+7.5*x*z^13
    7.5x*z^13 + 15.0x*y^2

    julia> @fastfunc "foo" p
    foo (generic function with 1 method)

    julia> foo(1.0,2.0,3.0)
    1.19574825e7
    ```
    """
macro fastfunc(functionname::String,polynomial)
	polystr = string(eval(polynomial))
	vars = eval(polynomial).x.vars
	# insert "*" in between numbers and variables
	regex = Regex("(\\d)" * replace(string(tuple(vars...)),", ","|"))
	polystr = replace(polystr,regex,s"\1*\2")
	# replace "+ -" with "- "
	polystr = replace(polystr,r"([+])(\s)([-])",s"\3 ")
	# create expression for function definition
	expr = parse(functionname*string(tuple(vars...))*" = @fastmath "*polystr)
	return eval(expr)
end
