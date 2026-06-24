# Limitations

RSDP v0.1 is intentionally narrow.

- It validates exact primal feasibility, not general optimality or infeasibility.
- It supports rational arithmetic, not algebraic number fields.
- Rational recovery is heuristic; certificate checking is exact.
- Real feasibility does not imply rational feasibility.
- Feasible sets on proper faces may require certified facial reduction.
- Exact PSD and affine algorithms are dense and intended for modest problems.
- Weighted SOS models can be validated after MOI bridging, but a readable polynomial
  identity requires metadata that may not be available from stable public APIs.
- Numerical solvers are replaceable hint generators and are never proof authorities.

When a feature is unsupported, RSDP should return a precise failure instead of silently
changing the model or weakening the claim.
