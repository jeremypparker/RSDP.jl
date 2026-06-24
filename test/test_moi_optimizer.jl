using Hypatia
using JuMP
using Test

const MOI = JuMP.MOI

struct OptimizerMissingOracle <: RSDP.AbstractNumericalOracle end

function _rsdp_optimizer_factory(; oracle = RSDP.HypatiaOracle())
    return () -> RSDP.Optimizer(
        oracle = oracle,
        max_denominator = big(10)^6,
    )
end

@testset "RSDP MOI optimizer" begin
    Q = Rational{BigInt}

    @testset "validated nonnegative feasibility" begin
        model = JuMP.GenericModel{Q}(_rsdp_optimizer_factory())
        @variable(model, x)
        @constraint(model, x == Q(1, 2))
        @constraint(model, x >= 0)
        @objective(model, Min, 2x + 1)

        optimize!(model)

        @test termination_status(model) == MOI.LOCALLY_SOLVED
        @test primal_status(model) == MOI.FEASIBLE_POINT
        @test result_count(model) == 1
        @test RSDP.validation_status(model) == RSDP.VALIDATED_PRIMAL_FEASIBLE
        @test RSDP.certificate(model) !== nothing
        @test RSDP.validation_report(model).ok
        @test RSDP.oracle_result(model).status ==
              RSDP.NUMERICAL_SOLVED_NOT_VALIDATED
        @test isempty(filter(contains("failed"), RSDP.diagnostics(model)))
        @test value(x) == Q(1, 2)
        @test JuMP.objective_value(model) == Q(2)
        @test raw_status(model) == "RSDP validated primal feasibility"
    end

    @testset "validated PSD triangle feasibility" begin
        model = JuMP.GenericModel{Q}(_rsdp_optimizer_factory())
        @variable(model, X[1:2, 1:2], Symmetric)
        @constraint(model, X in PSDCone())
        @constraint(model, X[1, 1] == Q(1))
        @constraint(model, X[1, 2] == Q(1, 2))
        @constraint(model, X[2, 2] == Q(1))

        optimize!(model)

        @test RSDP.validation_report(model).ok
        @test RSDP.certificate(model) !== nothing
        @test value(X[1, 1]) == Q(1)
        @test value(X[1, 2]) == Q(1, 2)
        @test value(X[2, 2]) == Q(1)
    end

    @testset "missing numerical implementation fails cleanly" begin
        model = JuMP.GenericModel{Q}(
            _rsdp_optimizer_factory(; oracle = OptimizerMissingOracle()),
        )
        @variable(model, x)
        @constraint(model, x == Q(1, 2))

        optimize!(model)

        @test termination_status(model) == MOI.NUMERICAL_ERROR
        @test primal_status(model) == MOI.NO_SOLUTION
        @test result_count(model) == 0
        @test RSDP.validation_status(model) == RSDP.NUMERICAL_ORACLE_FAILED
        @test isnothing(RSDP.certificate(model))
        @test any(contains("not loaded"), RSDP.diagnostics(model))
    end

    @testset "configuration attributes and strict scalar type" begin
        optimizer = RSDP.Optimizer()
        MOI.set(optimizer, RSDP.MaxDenominatorAttribute(), 1234)
        @test MOI.get(optimizer, RSDP.MaxDenominatorAttribute()) == 1234
        @test MOI.get(optimizer, RSDP.NumericalOracleAttribute()) isa
              RSDP.HypatiaOracle
        @test MOI.get(optimizer, RSDP.ExactificationPolicyAttribute()) isa
              RSDP.ErrorOnInexact
        @test !MOI.supports_constraint(
            optimizer,
            MOI.ScalarAffineFunction{Float64},
            MOI.EqualTo{Float64},
        )
    end
end
