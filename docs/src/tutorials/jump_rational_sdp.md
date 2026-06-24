# JuMP rational SDP

This example keeps every model coefficient rational. Hypatia is used only to
produce a floating-point hint; RSDP recovers a rational point and checks it
independently with exact arithmetic.

```julia
using Hypatia
using JuMP
using RSDP

const Q = Rational{BigInt}

model = JuMP.GenericModel{Q}()
@variable(model, X[1:3, 1:3], Symmetric)
@constraint(model, X in PSDCone())

target = Q[
    1    1//2 1//3
    1//2 1    1//4
    1//3 1//4 1
]
for column in 1:3, row in 1:column
    @constraint(model, X[row, column] == target[row, column])
end

extracted = RSDP.extract_moi(JuMP.backend(model))
result = RSDP.validate_with_oracle(
    extracted.problem,
    RSDP.HypatiaOracle();
    max_denominator=big(10)^6,
)

@assert result.certificate !== nothing
@assert result.report.ok
```

`result.report.ok` proves primal feasibility of the extracted rational conic
problem. It does not prove optimality: that would require an independently
checked exact dual certificate.

The complete tested source is
[`examples/jump_rational_sdp.jl`](https://github.com/jeremypparker/RSDP.jl/blob/master/examples/jump_rational_sdp.jl).

## Troubleshooting

If Hypatia solves but RSDP validation fails:

1. Increase `max_denominator`.
2. Try a less boundary-degenerate example.
3. Check that model coefficients are exact rationals; floating coefficients are
   rejected by default.
4. Inspect `result.diagnostics` and the rational-recovery diagnostics.
5. Consider whether certified facial reduction is needed.

Numerical success does not guarantee successful rational recovery.

