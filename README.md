
# SphericalHarmonicExpansions

| **Build Status** |
|:----------------:|
| [![CI](https://github.com/hofmannmartin/SphericalHarmonicExpansions.jl/workflows/CI/badge.svg)](https://github.com/hofmannmartin/SphericalHarmonicExpansions.jl/actions?query=workflow%3ACI) |
| [![codecov.io](http://codecov.io/gh/hofmannmartin/SphericalHarmonicExpansions.jl/coverage.svg?branch=master)](http://codecov.io/gh/hofmannmartin/SphericalHarmonicExpansions.jl?branch=master) |

The purpose of this package is to provide methods to numerically handle real spherical harmonics expansions in Cartesian coordinates.

## Table of Contents

- [1. Mathematical Background](#1-mathematical-background)
  - [1.1. Definition of the Spherical Harmonics](#11-definition-of-the-spherical-harmonics)
  - [1.2. Spherical Harmonics Expansions](#12-spherical-harmonics-expansions)
- [2. Usage](#2-usage)
  - [2.1. Polynomial Representation of the Spherical Harmonics](#21-polynomial-representation-of-the-spherical-harmonics)
  - [2.2. Polynomial Representation of the Spherical Harmonics Expansions](#22-polynomial-representation-of-the-spherical-harmonics-expansions)
  - [2.3. Transformation of Expansion Coefficients under Translation](#23-transformation-of-expansion-coefficients-under-translation)
  - [2.4. Transformation of Expansion Coefficients under Rotation](#24-transformation-of-expansion-coefficients-under-rotation)
  - [2.5. Numerical Evaluation](#25-numerical-evaluation)
- [3. Further Reading](#3-further-reading)

## 1. Mathematical Background

### 1.1. Definition of the Spherical Harmonics

The normalized real spherical harmonics on the unit sphere are defined by 
```math
Y_{l,m}(\vartheta,\varphi) := 
\begin{cases}
\sqrt{2}K_{l,m} \cos(m\varphi)P_{l,m}(\cos\vartheta) & m > 0\\
\sqrt{2}K_{l,m} \sin(-m\varphi)P_{l,-m}(\cos\vartheta) & m < 0\\
K_{l,m}P_{l,m}(\cos \vartheta) & m = 0
\end{cases},
``` 



where $`l\in\mathbb{N}_0`$, $`m\in [-l,l]`$, $`\theta`$ and $`\phi`$ are the spherical angular coordinates, 
```math
K_{l,m} = \sqrt{\frac{(2l+1)(l-|m|)!}{4\pi(l+|m|)!}},
```
is the normalization factor and
```math
P_{l,m}(x) = (1-x^2)^{\frac{m}{2}}\frac{d^m}{dx^m}\left(P_l(x)\right),
```

are the associated Legendre polynomials which can be derived from the Legendre polynomials
```math
P_l(x) = \frac{1}{2^ll!}\frac{d^l}{dx^l}\left[(x^2-1)^l\right].
```



Note that you will also find a convention in literature, where the  $`Y_{l,m}`$   are scaled by $`(-1)^m`$  . 

### 1.2. Spherical Harmonics Expansions
Each function  $f:\Omega \rightarrow \mathbb R$   satisfying Laplace's equation  $\Delta f = 0$   in a region  $\Omega\subseteq\mathbb R^3$   can be written as a spherical harmonic expansion
```math
f(\mathbf r) = \sum_{l=0}^{\infty}\sum_{m=-l}^l c_{l,m} r^l Y_l^m{\left(\frac{1}{r}\, \mathbf r\right)},
```

<div align="center"></div>

for all  $\mathbf a\in\Omega$  , where  $\mathbf c_{l,m}\in\mathbb R^3$   denote the spherical coefficients and  $r=\Vert \mathbf r \Vert_2$  . 

The term 
```math
r^l Y_l^m{\left(\frac{1}{r}\, \mathbf r\right)}
```

<div align="center"></div>

can be transformed from spherical to Cartesian coordinates, where it can be expressed as a homogeneous polynomial of degree  $l$  .

## 2. Usage
### 2.1. Polynomial Representation of the Spherical Harmonics
Generate a `MultivariatePolynomials.Polynomial` representation of
```math
Y_l^m{\left(\frac{1}{r}\, \mathbf r\right)}
```

<div align="center"></div>

in variables `α`, `β`, and `γ` on the unit sphere by

```julia
using SphericalHarmonics
@polyvar α β γ
l = 7 
m = -2

p = ylm(l,m,α,β,γ)
63.28217501963252αβγ⁵ - 48.67859616894809αβγ³ + 6.63799038667474αβγ
```

The polynomial representation of
```math
r^l Y_l^m{\left(\frac{1}{r}\, \mathbf r\right)}
```

<div align="center"></div>

in variables `x`, `y`, and `z` on  $\mathbb R^3$   can be obtained by

```julia
@polyvar x y z

p = rlylm(l,m,x,y,z)
6.63799038667474x⁵yz + 13.27598077334948x³y³z - 35.40261539559861x³yz³ + 6.63799038667474xy⁵z - 35.40261539559861xy³z³ + 21.24156923735917xyz⁵
```

### 2.2. Polynomial Representation of the Spherical Harmonics Expansions
In case where a function is equal to or can be approximated by a **finite** Spherical harmonic expansion
```math
\sum_{l=0}^{L}\sum_{m=-l}^l c_{l,m} r^l Y_l^m{\left(\frac{1}{r}\, \mathbf r\right)},
```

<div align="center"></div>

with  $L \in \mathbb N$   its multivariate polynomial representation has finite degree.

Coefficents  $c_{l,m}$   can be initialized and populated by `c[l,m] = 42.0`.

```julia
L = 2
c = SphericalHarmonicCoefficients(L)
c[0,0] = 42.0 #c₀₀
c[2,-1] = -1.0 #c₂₋₁
c[2,1] = 2.0 #c₂₁
```
Internally, the coefficients are lexicographically stored in a vector (`c[0,0]`, `c[1,-1]`, `c[1,0]`, `c[1,1]`, `c[2,-2]`, ...). So the above initialization is equivalent to
```julia
C = [42.0,0,0,0,0,-1,0,2,0]
c = SphericalHarmonicCoefficients(C)
```
Note that `SphericalHarmonicCoefficients(C)` will throw an error if `length(C)` is not  $(L+1)^2$   for some  $L\in\mathbb{N}$. From there on the corresponding polynomial  representation in Cartesian coordinates `x`, `y`, and `z` can be obtained by 
```julia
@polyvar x y z

f = sphericalHarmonicsExpansion(c,x,y,z)
11.847981254502882 - 1.0925484305920792yz + 2.1850968611841584xz
```
Currently, expansions up to $L=66$ are supported.

### 2.3. Transformation of Expansion Coefficients under Translation

If we change from a coordinate system with coordinates `x`, `y`, and `z` into a translated one with new coordinates `u = x + tx`, `v = y + ty`, and `w = z + tz` we need transformed coefficients to express the expansion in these new coordinates. To this end, we can do 

```julia
@polyvar u v w
translationVector = [0,0,1.0] # [tx,ty,tz]

cTranslated = translation(c,translationVector)
sphericalHarmonicsExpansion(cTranslated,u,v,w)
11.847981254502878 - 1.0925484305920792v + 2.1850968611841584u - 1.092548430592079vw + 2.185096861184158uw
```

### 2.4. Transformation of Expansion Coefficients under Rotation

Furthermore, any rotation of the coordinate system can be expressed with rotated coefficients. Rotation of a coordinate system with coordinates `x`, `y`, and `z` into a rotated system with coordinates `(r₁, r₂, r₃) = R * [x, y, z]`, using a rotation matrix `R = RotZYZ(eulerAngles...)` provided by [Rotations.jl](https://github.com/JuliaGeometry/Rotations.jl), can be obtained by a rotation of the coefficients with

```julia
@polyvar r₁ r₂ r₃
rotationAngles = [π/2, -π, 0] # ZYZ-convention (Euler angles)

cRotated = rotation(c, rotationAngles)
sphericalHarmonicsExpansion(cRotated, r₁,r₂,r₃)
11.847981254502882 + 2.185096861184156r₂r₃ - 1.0925484305920796r₁r₃
```

Combining the rotation with a point reflection and translation, any reflections of the coordinate system can be represented. For example, a reflection of the coordinate system with coordinates `x`, `y`, and `z` in the xz-plane yielding a coordinates `x`, `y⁻ = -y`, and `z` can be expressed in the coefficients using

```julia
@polyvar y⁻
        
cReflected = rotation(pointReflection(c), [0,pi,0]);
sphericalHarmonicsExpansion(cReflected, x, y⁻, z)
11.847981254502882 + 1.0925484305920798y⁻z + 2.185096861184156xz
```

### 2.5. Numerical Evaluation

If you want to evaluate  $f$   at a specific point you can use the standard interface of `MultivariatePolynomials`

```julia
f(x=>0.5, y=>-1.0, z=>0.25)
12.394255469798921
f((x,y,z)=>(0.5,-1.0,0.25))
12.394255469798921
```

In case where you want to evaluate  $f$   for a large number of points you might run into performance issues. To this end we provide two methods to dynamically generate fast evaluating functions. Either use

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

## 3. Further Reading

For more informations on the `MultivariatePolynomials` package please visit the project page on [github](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl).
