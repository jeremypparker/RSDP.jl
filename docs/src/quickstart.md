# Quick start

Create a problem in standard conic form:

```math
A x = b,\qquad x \in K.
```

```julia
using RSDP

Q = Rational{BigInt}
A = Q[1 1 0 0 0]
b = Q[1]
cones = AbstractConeBlock[
    NonnegativeConeBlock(2),
    PSDTriangleConeBlock(2),
]
problem = ExactConicProblem(A, b, cones)

x = Q[1//2, 1//2, 1, 0, 1]
certificate = make_primal_certificate(problem, x)
report = check_certificate(problem, certificate)
@assert report.ok
```

This proves only that the supplied point satisfies the rational affine equations and
belongs to the product cone.

For an approximate point, use `recover_primal_certificate`. Recovery can fail even
when a nearby real feasible point exists; failure is reported rather than promoted to
a mathematical claim.

With the optional Hypatia dependency loaded, `validate_with_oracle(problem,
HypatiaOracle())` performs numerical hint generation, rational recovery, and
the same independent exact certificate check. It validates primal feasibility,
not optimality.
