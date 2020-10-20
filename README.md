# SphericalHarmonics

| **Build Status** |
|:----------------:|
| [![Build Status](https://travis-ci.org/hofmannmartin/SphericalHarmonics.jl.svg?branch=master)](https://travis-ci.org/hofmannmartin/SphericalHarmonics.jl) |
| [![codecov](https://codecov.io/gh/hofmannmartin/SphericalHarmonics.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/hofmannmartin/SphericalHarmonics.jl) |

The purpose of this package is to provide methods to numerically handle real spherical harmonics expansions in Cartesian coordinates.

## Table of Contents

- [Mathematical Background](#mathematical-background)
  - [Definition of the Spherical Harmonics](#definition-of-the-spherical-harmonics)
  - [Spherical Harmonics Expansions](#spherical-harmonics-expansions)
- [Usage](#usage)
  - [Polynomial Representation of the Spherical Harmonics](#polynomial-representation-of-the-spherical-harmonics)
  - [Polynomial Representation of the Spherical Harmonics Expansions](#polynomial-representation-of-the-spherical-harmonics-expansions)
  - [Transformation of Expansion Coefficients under Translation](#transformation-of-expansion-coefficients-under-translation)
  - [Numerical Evaluation](#numerical-valuation)
- [Further Reading](#further-reading)

## Mathematical Background

### Definition of the Spherical Harmonics

The normalized real spherical harmonics on the unit sphere are defined by 
<!-- $$
Y_{l,m}(\vartheta,\varphi) := 
\begin{cases}
\sqrt{2}K_{l,m} \cos(m\varphi)P_{l,m}(\cos\vartheta) & m > 0\\
\sqrt{2}K_{l,m} \sin(-m\varphi)P_{l,-m}(\cos\vartheta) & m < 0\\
K_{l,m}P_{l,m}(\cos \vartheta) & m = 0,
\end{cases}
$$ --> 

<div align="center"><img src="https://render.githubusercontent.com/render/math?math=Y_%7Bl%2Cm%7D(%5Cvartheta%2C%5Cvarphi)%20%3A%3D%20%0A%5Cbegin%7Bcases%7D%0A%5Csqrt%7B2%7DK_%7Bl%2Cm%7D%20%5Ccos(m%5Cvarphi)P_%7Bl%2Cm%7D(%5Ccos%5Cvartheta)%20%26%20m%20%3E%200%5C%5C%0A%5Csqrt%7B2%7DK_%7Bl%2Cm%7D%20%5Csin(-m%5Cvarphi)P_%7Bl%2C-m%7D(%5Ccos%5Cvartheta)%20%26%20m%20%3C%200%5C%5C%0AK_%7Bl%2Cm%7DP_%7Bl%2Cm%7D(%5Ccos%20%5Cvartheta)%20%26%20m%20%3D%200%2C%0A%5Cend%7Bcases%7D"></div>

where <!-- $l\in\mathbb{N}_0$ --> <img src="https://render.githubusercontent.com/render/math?math=l%5Cin%5Cmathbb%7BN%7D_0"> and <!-- $m\in [-l,l]$ --> <img src="https://render.githubusercontent.com/render/math?math=m%5Cin%20%5B-l%2Cl%5D">, <!-- $\theta$ --> <img src="https://render.githubusercontent.com/render/math?math=%5Ctheta"> and <!-- $\phi$ --> <img src="https://render.githubusercontent.com/render/math?math=%5Cphi"> are the spherical angular coordinates, 
<!-- $$
K_{l,m} = \sqrt{\frac{(2l+1)(l-|m|)!}{4\pi(l+|m|)!}},
$$ --> 

<div align="center"><img src="https://render.githubusercontent.com/render/math?math=K_%7Bl%2Cm%7D%20%3D%20%5Csqrt%7B%5Cfrac%7B(2l%2B1)(l-%7Cm%7C)!%7D%7B4%5Cpi(l%2B%7Cm%7C)!%7D%7D%2C"></div>

is the normalization factor and
<!-- $$
P_{l,m}(x) = (1-x^2)^{\frac{m}{2}}\frac{d^m}{dx^m}\left(P_l(x)\right),
$$ --> 

<div align="center"><img src="https://render.githubusercontent.com/render/math?math=P_%7Bl%2Cm%7D(x)%20%3D%20(1-x%5E2)%5E%7B%5Cfrac%7Bm%7D%7B2%7D%7D%5Cfrac%7Bd%5Em%7D%7Bdx%5Em%7D%5Cleft(P_l(x)%5Cright)%2C"></div>

are the associated Legendre polynomials which can be derived from the Legendre polynomials
<!-- $$
P_l(x) = \frac{1}{2^ll!}\frac{d^l}{dx^l}\left[(x^2-1)^l\right].
$$ --> 

<div align="center"><img src="https://render.githubusercontent.com/render/math?math=P_l(x)%20%3D%20%5Cfrac%7B1%7D%7B2%5Ell!%7D%5Cfrac%7Bd%5El%7D%7Bdx%5El%7D%5Cleft%5B(x%5E2-1)%5El%5Cright%5D."></div>

Note that you will also find a convention in literature, where the <!-- $Y_{l,m}$ --> <img src="https://render.githubusercontent.com/render/math?math=Y_%7Bl%2Cm%7D"> are scaled by <!-- $(-1)^m$ --> <img src="https://render.githubusercontent.com/render/math?math=(-1)%5Em">. 

### Spherical Harmonics Expansions
Each function <!-- $f:\Omega \rightarrow \mathbb R$ --> <img src="https://render.githubusercontent.com/render/math?math=f%3A%5COmega%20%5Crightarrow%20%5Cmathbb%20R"> satisfying Laplace's equation <!-- $\Delta f = 0$ --> <img src="https://render.githubusercontent.com/render/math?math=%5CDelta%20f%20%3D%200"> in a region <!-- $\Omega\subseteq\mathbb R^3$ --> <img src="https://render.githubusercontent.com/render/math?math=%5COmega%5Csubseteq%5Cmathbb%20R%5E3"> can be written as a spherical harmonic expansion
<!-- $$
f(\mathbf r) = \sum_{l=0}^{\infty}\sum_{m=-l}^l c_{l,m} r^l Y_l^m{\left(\frac{1}{r}\, \mathbf r\right)},
$$ --> 

<div align="center"><img src="https://render.githubusercontent.com/render/math?math=f(%5Cmathbf%20r)%20%3D%20%5Csum_%7Bl%3D0%7D%5E%7B%5Cinfty%7D%5Csum_%7Bm%3D-l%7D%5El%20c_%7Bl%2Cm%7D%20r%5El%20Y_l%5Em%7B%5Cleft(%5Cfrac%7B1%7D%7Br%7D%5C%2C%20%5Cmathbf%20r%5Cright)%7D%2C"></div>

for all <!-- $\mathbf a\in\Omega$ --> <img src="https://render.githubusercontent.com/render/math?math=%5Cmathbf%20a%5Cin%5COmega">, where <!-- $\mathbf c_{l,m}\in\mathbb R^3$ --> <img src="https://render.githubusercontent.com/render/math?math=%5Cmathbf%20c_%7Bl%2Cm%7D%5Cin%5Cmathbb%20R%5E3"> denote the spherical coefficients and <!-- $r=\Vert \mathbf r \Vert_2$ --> <img src="https://render.githubusercontent.com/render/math?math=r%3D%5CVert%20%5Cmathbf%20r%20%5CVert_2">. 

The term 
<!-- $$
r^l Y_l^m{\left(\frac{1}{r}\, \mathbf r\right)}
$$ --> 

<div align="center"><img src="https://render.githubusercontent.com/render/math?math=r%5El%20Y_l%5Em%7B%5Cleft(%5Cfrac%7B1%7D%7Br%7D%5C%2C%20%5Cmathbf%20r%5Cright)%7D"></div>

can be transformed from from spherical to Cartesian coordinates, where is can be expressed as a homogeneous polynomial of degree <!-- $l$ --> <img src="https://render.githubusercontent.com/render/math?math=l">.

## Usage
### Polynomial Representation of the Spherical Harmonics
Generate a `MultivariatePolynomials.Polynomial` representation of
<!-- $$
Y_l^m{\left(\frac{1}{r}\, \mathbf r\right)}
$$ --> 

<div align="center"><img src="https://render.githubusercontent.com/render/math?math=Y_l%5Em%7B%5Cleft(%5Cfrac%7B1%7D%7Br%7D%5C%2C%20%5Cmathbf%20r%5Cright)%7D"></div>

in variables `x̂`, `ŷ`, and `ẑ` on the unit sphere by

```julia
using SphericalHarmonics
@polyvar x̂ ŷ ẑ
l = 7 
m = -2

p = ylm(l,m,x̂,ŷ,ẑ)
63.28217501963252x̂ŷẑ⁵ - 48.67859616894809x̂ŷẑ³ + 6.63799038667474x̂ŷẑ
```

The polynomial representation of
<!-- $$
r^l Y_l^m{\left(\frac{1}{r}\, \mathbf r\right)}
$$ --> 

<div align="center"><img src="https://render.githubusercontent.com/render/math?math=r%5El%20Y_l%5Em%7B%5Cleft(%5Cfrac%7B1%7D%7Br%7D%5C%2C%20%5Cmathbf%20r%5Cright)%7D"></div>

in variables `x`, `y`, and `z` on <!-- $\mathbb R^3$ --> <img src="https://render.githubusercontent.com/render/math?math=%5Cmathbb%20R%5E3"> can be obtained by

```julia
@polyvar x y z

p = rlylm(l,m,x,y,z)
6.63799038667474x⁵yz + 13.27598077334948x³y³z - 35.40261539559861x³yz³ + 6.63799038667474xy⁵z - 35.40261539559861xy³z³ + 21.24156923735917xyz⁵
```

### Polynomial Representation of the Spherical Harmonics Expansions
In case where a function is equal to or can be approximated by a **finite** Spherical harmonic expansion
<!-- $$
\sum_{l=0}^{L}\sum_{m=-l}^l c_{l,m} r^l Y_l^m{\left(\frac{1}{r}\, \mathbf r\right)},
$$ --> 

<div align="center"><img src="https://render.githubusercontent.com/render/math?math=%5Csum_%7Bl%3D0%7D%5E%7BL%7D%5Csum_%7Bm%3D-l%7D%5El%20c_%7Bl%2Cm%7D%20r%5El%20Y_l%5Em%7B%5Cleft(%5Cfrac%7B1%7D%7Br%7D%5C%2C%20%5Cmathbf%20r%5Cright)%7D%2C"></div>

with <!-- $L \in \mathbb N$ --> <img src="https://render.githubusercontent.com/render/math?math=L%20%5Cin%20%5Cmathbb%20N"> its multivariate polynomial representation has finite degree.

Coefficents <!-- $c_{l,m}$ --> <img src="https://render.githubusercontent.com/render/math?math=c_%7Bl%2Cm%7D"> can be initialized and populated by by `c[l,m] = 42.0`.

```julia
L = 2
c = SphericalHarmonicCoefficients(L)
c[0,0] = 42.0 #c₀₀
c[2,-1] = -1.0 #c₂₋₁
c[2,1] = 2.0 #c₂₁
```
Internally the coefficients are lexicographically stored in a vector (`c[0,0]`, `c[1,-1]`, `c[1,0]`, `c[1,1]`, `c[2,-2]`, ...). So the above initialization is equivalent to
```julia
C = [42.0,0,0,0,0,-1,0,2,0]
c = SphericalHarmonicCoefficients(C)
f = sphericalHarmonicsExpansion(c,x,y,z)
2.1850968611841584xz + -1.0925484305920792yz + 11.847981254502882
```
Note that `SphericalHarmonicCoefficients(C)` will throw an error if `length(C)` is not <!-- $(L+1)^2$ --> <img src="https://render.githubusercontent.com/render/math?math=(L%2B1)%5E2"> for some <!-- $L\in\mathbb{N}$ --> <img src="https://render.githubusercontent.com/render/math?math=L%5Cin%5Cmathbb%7BN%7D">. From there on the corresponding polynomial  representation in cartesian coordinates `x`, `y`, and `z` can be obtained by 
```julia
@polyvar x y z

f = sphericalHarmonicsExpansion(c,x,y,z)
2.1850968611841584xz - 1.0925484305920792yz + 11.847981254502882
```

### Transformation of Expansion Coefficients under Translation

If we change from a coordinate sytsem with coordinates `x`, `y`, and `z` into a translated one with new coordinates `u = x + tx`, `v = y + ty`, and `w = z + tz` we need transformed coefficients to express the expansiion in these new coordinates. To this end we can do 

```julia
@polyvar u v w
translationVector = [0,0,1.0] # [tx,ty,tz]

cTranslated = translation(c,translationVector)
sphericalHarmonicsExpansion(cTranslated,u,v,w)
2.1850968611841584uw - 1.0925484305920792vw + 2.1850968611841584u - 1.0925484305920792v + 11.847981254502878
```

### Numerical Evaluation

If you want to evaluate <!-- $f$ --> <img src="https://render.githubusercontent.com/render/math?math=f"> at a specific point you can use the standard interface of `MultivariatePolynomials`

```julia
f(x=>0.5, y=>-1.0, z=>0.25)
12.394255469798921
f((x,y,z)=>(0.5,-1.0,0.25))
12.394255469798921
```

In case where you want to evaluate <!-- $f$ --> <img src="https://render.githubusercontent.com/render/math?math=f"> for a large number of points you might run into performance issues. To this end we provide two methods to dynamically generate fast evaluating functions. Either use

```julia
g = @fastfunc f
g(0.5,-1.0,0.25)
12.394255469798921
```
which has moderate generation overhead. Usage from within local scope requires `Base.invokelatest(foo, 1.0,2.0,3.0)` instead of `foo(1.0,2.0,3.0)` to avoid issue #4. Or use

```julia
h = fastfunc(f)
h(0.5,-1.0,0.25)
12.394255469798921
```
which uses `GeneralizedGenerated` for function generation and comes with a significant overhead.

## Further Reading

For more informations on the `MultivariatePolynomials` package please visit the project page on [github](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl).