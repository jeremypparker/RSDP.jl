# RSDP Optimizer

`RSDP.Optimizer` provides a direct JuMP workflow for validated primal
feasibility:

```julia
using Hypatia
using JuMP
using RSDP

const Q = Rational{BigInt}

model = JuMP.GenericModel{Q}(() -> RSDP.Optimizer(
    oracle = RSDP.HypatiaOracle(),
    max_denominator = big(10)^6,
))

@variable(model, X[1:2, 1:2], Symmetric)
@constraint(model, X in PSDCone())
@constraint(model, X[1, 1] == Q(1))
@constraint(model, X[1, 2] == Q(1, 2))
@constraint(model, X[2, 2] == Q(1))

optimize!(model)

@assert RSDP.validation_status(model) == RSDP.VALIDATED_PRIMAL_FEASIBLE
@assert RSDP.validation_report(model).ok
cert = RSDP.certificate(model)
@assert cert !== nothing
```

The optimizer wraps the existing pipeline:

```text
MOI model -> exact extraction -> numerical oracle -> rational recovery
          -> independent exact certificate check
```

This validates primal feasibility only. `RSDP.Optimizer` deliberately reports
`MOI.LOCALLY_SOLVED`, not `MOI.OPTIMAL`, because no exact dual optimality
certificate is produced.

Hypatia generates only a floating-point hint. The final certificate and
validation report come from RSDP's rational recovery and exact checker, which
are independent of Hypatia. Values returned by `value` are reconstructed from
the exact certificate rather than copied from Hypatia's raw solution.

Use exact model data, normally with
`JuMP.GenericModel{Rational{BigInt}}`. Floating-point functions and sets are not
advertised as supported, and extraction rejects floating coefficients by
default. An explicit exactification policy is required to opt into conversion.

Inspect `RSDP.diagnostics(model)` if extraction, numerical solution, or rational
recovery fails.
