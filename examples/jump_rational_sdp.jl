module JuMPRationalSDPExample

using Hypatia
using JuMP
using RSDP

const Q = Rational{BigInt}

function run()
    model = JuMP.GenericModel{Q}()
    @variable(model, X[1:3, 1:3], Symmetric)
    @constraint(model, X in PSDCone())

    target = Q[
        1 1//2 1//3
        1//2 1 1//4
        1//3 1//4 1
    ]
    for column in 1:3, row in 1:column
        @constraint(model, X[row, column] == target[row, column])
    end

    extracted = RSDP.extract_moi(JuMP.backend(model))
    result = RSDP.validate_with_oracle(
        extracted.problem,
        RSDP.HypatiaOracle();
        max_denominator = big(10)^6,
    )
    @assert result.certificate !== nothing
    @assert result.report.ok
    return result
end

const result = run()

end

JuMPRationalSDPExample.result
