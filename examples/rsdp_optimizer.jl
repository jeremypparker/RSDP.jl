module RSDPOptimizerExample

using Hypatia
using JuMP
using RSDP

const Q = Rational{BigInt}

function run()
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
    @assert RSDP.certificate(model) !== nothing
    return RSDP.validation_report(model)
end

const report = run()

end

RSDPOptimizerExample.report
