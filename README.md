# FunQuad.jl

This package enables globally h-adaptive quadrature for smooth functions and
functions with endpoint singularities. The usage is similar to QuadGK.jl,
however, instead of passing an interval with breakpoints, the caller will pass
multiple function spaces, each with an interval and the kinds of functions to
expect on that interval.

For example, the behavior of QuadGK.jl can be reproduced as follows
```julia
using FunQuad
funquadh(x -> sin(10x), FunSpace(0, 1))
```

The default functions for a `FunSpace` are `GaussKronrod`, which are polynomials
orthogonal w.r.t. the Legendre weight function.
To specify other kinds of functions, pass a function kind to the first argument
of `FunSpace`. For example, integrating logarithmic singularities as follows
```julia
funquadh(x -> log(x) *sin(8x) + x^3, FunSpace(LogAtMin(), 0, 1))
```

The caller should add breakpoints to the initial domain so that there is at most
one kind of singularity in any interval. For example (not tested)
```julia
funquadh(x -> log(x) *sin(8x) + x^3*log(1-x), FunSpace(LogAtMin(), 0, 0.5), FunSpace(LogAtMax(), 0.5, 1))
```

This package will attempt to evaluate the integrand in the interior of the
interval, not at the endpoints. If available, it will use Kronrod rules to
compute embedded error estimates (otherwise refines the panel). The singularity
should be included in the integrand, not separated into a integration measure.
The quadrature rules attempt to integrate functions of the form smooth plus
singular times smooth to attain faster convergence for a wide class of
functions.

Currently, the following function spaces are available:
- `GaussKronrod()`: integrates smooth functions with embedded error estimates.
 (Default function space)
- `GaussLegendre()`: integrates smooth functions with refined error estimates.
- `LogAtMin()`: integrates logarithmic singularities at left endpoints. Refined
  error estimates
- `LogAtMax()`: integrates logarithmic singularities at right endpoints. Refined
  error estimates

TODO: add a `FunSpace` for changes of variables, such as
`FunSpace(Transformed(GaussKronrod, x -> sqrt(x), y -> y^2, y -> 2y), 0, 1)`

TODO: add power laws with GeneralizedChebyshevQuadrature.jl

TODO: add Jacobi weight functions with QuadJGK.jl