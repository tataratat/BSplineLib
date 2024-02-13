# BSplineLib
BSpline and NURBS library written in modern c++ for seamless integration to
[splinepy](https://github.com/tataratat/splinepy).
This is a forked project of [SplineLib](https://github.com/SplineLib/SplineLib).

## Why a fork?
`SplineLib` implements [this paper](https://doi.org/10.1016/j.advengsoft.2020.102826) under MIT License, where they take "basis function as a data structure" approach; each basis function is an object and holds precomputed denominators of Cox-de Boor recursion formula.
Such basis functions are created in fully recursive manner and would be also retrieved in fully recursive manner during evaluations.
While this works fine with low-order splines, it posed some practical issues with high-order splines.
Increasing orders/degrees results in increased depth of recursion, which caused high memory consumption and poor performance.
At the same time, bookkeeping of the basis functions became cumbersome in dynamic situations we aim to support in `splinepy`.
Due to aforementioned inconveniences in our applications, this project opted for `A2.2` & `A2.3` algorithms from `The NURBS Book`.

## What can `BSplineLib` do / does differently?
- Basis function evaluations using `The NURBS BOOK - A2.2`
- Basis function derivative evaluations using `The NURBS BOOK - A2.3`
- Direct contiguous buffer access for evaluation input and output
- Less use of `NamedType`
- Control points (VectorSpace) creation without copy, based on a contiguous memory view
- Efficient tensor product using adapted version of bezman's [RecursiveCombine](https://github.com/tataratat/bezman/blob/76ff2d9d6f29ea333b5f727a0ecee3d43e07e6e5/src/utils/algorithms/recursive_combine.hpp#L39)
- Tests implemented in `splinepy`
