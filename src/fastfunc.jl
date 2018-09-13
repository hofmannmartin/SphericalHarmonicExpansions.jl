"""
    @fastfunc functionname[::String] polynomial[::Polynomial{true}]

    Defines a function with name `functionname` for `polynomial` for fast numerical evaluation.
    Note that the polynomial is transformed using `@fastmath`, which calls functions that
    may violate strict IEEE semantics.

    # Examples
    ```
    julia> using MultivariatePolynomials

    julia> p = 15.0*x*y^2+7.5*x*z^13;

    julia> @fastfunc "foo" p;

    julia> foo(1.0,2.0,3.0)
    1.19574825e7
    ```
    """
macro fastfunc(polynomial)
	return quote
		local polystr = string($(esc(polynomial)))
		local vars = string(tuple(variables($(esc(polynomial)))...))
		# insert "*" in between numbers and variables
		local regex = Regex("(\\d)" * replace(vars,", " => "|"))
		polystr = replace(polystr,regex => s"\1*\2")
		# insert "*" in between variables "xy" -> "x*y"
		local regex_xy = Regex(replace(vars,", " => "|")*replace(vars,", " => "|"))
		for i=1:length(variables($(esc(polynomial))))-1
			polystr = replace(polystr,regex_xy => s"\1*\2")
		end
		# replace "+ -" with "- "
		polystr = replace(polystr,r"([+])(\s)([-])" => s"\3 ")
		# create expression for function definition

		eval(Meta.parse(vars*" -> @fastmath "*polystr))
	end
end
