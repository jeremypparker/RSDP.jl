# Ordinary sum of squares

Consider the strictly SOS polynomial

```math
p(x) = 1 + x^2 + x^4.
```

SumOfSquares currently constructs a `Float64` JuMP model, so this tutorial does
not silently treat its floating bridge data as exact. It solves that model with
Hypatia, then reconstructs the simple rational Gram formulation using only
public APIs:

```math
p(x) = [1,x,x^2]^\mathsf{T} I [1,x,x^2].
```

```julia
using DynamicPolynomials
using Hypatia
using JuMP
using LinearAlgebra: I
using RSDP
using SumOfSquares

const Q = Rational{BigInt}
const MOI = JuMP.MOI

@polyvar x
p = 1 + x^2 + x^4

model = SOSModel(Hypatia.Optimizer)
@constraint(model, p in SOSCone())
set_silent(model)
optimize!(model)
@assert termination_status(model) in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL)

# Exact Gram SDP for p = [1, x, x^2]' * I * [1, x, x^2].
gram_model = JuMP.GenericModel{Q}()
@variable(gram_model, G[1:3, 1:3], Symmetric)
@constraint(gram_model, G in PSDCone())
for column in 1:3, row in 1:column
    @constraint(gram_model, G[row, column] == Q(row == column))
end

extracted = RSDP.extract_moi(JuMP.backend(gram_model))
result = RSDP.validate_with_oracle(
    extracted.problem,
    RSDP.HypatiaOracle();
    max_denominator=big(10)^6,
)

@assert result.certificate !== nothing
@assert result.report.ok
```

The first Hypatia solve is numerical. The second Hypatia call, inside
`validate_with_oracle`, is also only a hint generator. The final claim comes
from exact rational recovery, exact PSD checking, and the independent
certificate checker.

This example validates the exact Gram SDP and the displayed polynomial
identity. General extraction of a human-readable SOS identity from
SumOfSquares bridge metadata is not yet implemented, and private
SumOfSquares internals are deliberately not used.

Only primal feasibility is validated, not optimality. Floating coefficients
are rejected by default, and rational recovery may fail even after a successful
numerical solve.

The complete tested source is
[`examples/ordinary_sos.jl`](https://github.com/jeremypparker/RSDP.jl/blob/master/examples/ordinary_sos.jl).

## Troubleshooting

If Hypatia solves but RSDP validation fails:

1. Increase `max_denominator`.
2. Try a less boundary-degenerate example.
3. Check that the exact Gram model uses integers or rationals.
4. Inspect `result.diagnostics`.
5. Consider whether certified facial reduction is needed.
