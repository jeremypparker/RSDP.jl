module OrdinarySOSExample

using DynamicPolynomials
using Hypatia
using JuMP
using LinearAlgebra: I
using RSDP
using SumOfSquares

const Q = Rational{BigInt}
const MOI = JuMP.MOI

function run()
    @polyvar x
    polynomial = 1 + x^2 + x^4

    # SumOfSquares currently builds Float64 JuMP models. This solve is a
    # numerical existence check only, not the exact validation step.
    sos_model = SOSModel(Hypatia.Optimizer)
    @constraint(sos_model, polynomial in SOSCone())
    set_silent(sos_model)
    optimize!(sos_model)
    termination_status(sos_model) in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL) ||
        error("Hypatia did not solve the SumOfSquares model")

    # Reconstruct the exact public Gram formulation for
    # polynomial = [1, x, x^2]' * I * [1, x, x^2].
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
        max_denominator = big(10)^6,
    )
    @assert result.certificate !== nothing
    @assert result.report.ok

    basis = [one(x), x, x^2]
    exact_gram = Matrix{Q}(I, 3, 3)
    reconstructed = sum(
        exact_gram[row, column] * basis[row] * basis[column] for row in 1:3, column in 1:3
    )
    @assert reconstructed == polynomial
    return result
end

const result = run()

end

OrdinarySOSExample.result
