# QuasiPeriodicHelmholtz.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mjp98.github.io/QuasiPeriodicHelmholtz.jl/dev)
[![Build Status](https://github.com/mjp98/QuasiPeriodicHelmholtz.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mjp98/QuasiPeriodicHelmholtz.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mjp98/QuasiPeriodicHelmholtz.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mjp98/QuasiPeriodicHelmholtz.jl)

This is an experimental package to compute quasi-periodic Green's functions for the Helmholtz equation in 2D.

It implements the tail end summation correction to compute the quasi-periodic Green's function

```julia
    G = GreensFunction()
```

Dispatch is used to return derivatives and decompositions of a Green's function in the following way:

```julia
    G(x::SVector{2},y::SVector{2})
```

computes 

$$ G(x,y) $$


```julia
    G(x::SVector{2},y::SVector{2},dx::SVector{2})
```

computes 

$$ \left(dx \cdot \frac{\partial }{\partial x}\right) G(x,y) $$

```julia
    G(x::SVector{2},y::SVector{2},dx::SVector{2},dy::SVector{2})
```

computes 

$$ \left(dy \cdot \frac{\partial }{\partial y}\right) \left(dx \cdot \frac{\partial }{\partial x}\right) G(x,y) $$

Partial support is provided for a decomposition of the form

```julia
    G(x::SVector{2},y::SVector{2}) =\Sigma_{n \in \{0,1,2\}} G(n,x,y) (y-x)^(-n)/\pi + G(:log,x,y)\log(|y-x|)/\pi
```


where `G(n,x,y)` is analytic at `x=y`
