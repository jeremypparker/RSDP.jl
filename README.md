# RSDP.jl

[![Build Status](https://github.com/jeremypparker/RSDP.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/jeremypparker/RSDP.jl/actions/workflows/CI.yml)
[![Coverage](https://codecov.io/gh/jeremypparker/RSDP.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jeremypparker/RSDP.jl)

RSDP validates exact rational certificates for modest semidefinite and conic
feasibility problems. Numerical points are hints, never proofs.

The initial commit focuses on:

- rational affine equalities;
- zero, nonnegative, and positive-semidefinite triangle cone blocks;
- exact primal certificates and objective evaluation;
- rational recovery in an exact affine space;
- strict MathOptInterface extraction; and
- clear diagnostics for unsupported, inexact, and boundary cases.

```julia
using RSDP

problem = ExactConicProblem(A, b, cones; objective=c)
certificate = make_primal_certificate(problem, rational_point)
report = check_certificate(problem, certificate)
@assert report.ok
```

RSDP does **not** infer exact model data from floating-point literals by default. It
does not claim infeasibility or optimality without the corresponding exact dual
certificate. Automatic facial reduction, number fields, and human-readable weighted
SOS reconstruction are planned work.

See the [development documentation](https://jeremypparker.github.io/RSDP.jl/dev/) and
[`docs/src/manual/limitations.md`](docs/src/manual/limitations.md).
