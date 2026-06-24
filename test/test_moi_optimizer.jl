using Hypatia
using JuMP
using Test

const MOI = JuMP.MOI

struct OptimizerMissingOracle <: RSDP.AbstractNumericalOracle end
struct OptimizerRecoveryFailureOracle <: RSDP.AbstractNumericalOracle end
struct OptimizerUnsupportedConeOracle <: RSDP.AbstractNumericalOracle end

function RSDP.solve_oracle(
    ::RSDP.ExactConicProblem,
    ::OptimizerRecoveryFailureOracle;
    kwargs...,
)
    return RSDP.NumericalOracleResult(
        RSDP.NUMERICAL_SOLVED_NOT_VALIDATED;
        primal = [-1.0, 0.0],
    )
end

function RSDP.solve_oracle(
    ::RSDP.ExactConicProblem,
    ::OptimizerUnsupportedConeOracle;
    kwargs...,
)
    return RSDP.NumericalOracleResult(
        RSDP.UNSUPPORTED_CONE;
        diagnostics = ["unsupported test cone"],
    )
end

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
        @test MOI.get(optimizer, MOI.TerminationStatus()) == MOI.OPTIMIZE_NOT_CALLED
        @test MOI.get(optimizer, MOI.DualStatus()) == MOI.NO_SOLUTION
        @test MOI.get(optimizer, MOI.SolverName()) ==
              "RSDP validated primal feasibility"
        @test isnothing(MOI.get(optimizer, RSDP.CertificateAttribute()))
        @test isnothing(MOI.get(optimizer, RSDP.ValidationReportAttribute()))
        @test isnothing(MOI.get(optimizer, RSDP.OracleResultAttribute()))
        @test isempty(MOI.get(optimizer, RSDP.DiagnosticsAttribute()))

        @test MOI.supports(optimizer, MOI.Silent())
        MOI.set(optimizer, MOI.Silent(), false)
        @test !MOI.get(optimizer, MOI.Silent())

        MOI.set(optimizer, RSDP.MaxDenominatorAttribute(), 1234)
        @test MOI.get(optimizer, RSDP.MaxDenominatorAttribute()) == 1234
        @test_throws ArgumentError MOI.set(
            optimizer,
            RSDP.MaxDenominatorAttribute(),
            0,
        )
        @test_throws ArgumentError RSDP.Optimizer(; max_denominator = 0)

        missing_oracle = OptimizerMissingOracle()
        MOI.set(optimizer, RSDP.NumericalOracleAttribute(), missing_oracle)
        @test MOI.get(optimizer, RSDP.NumericalOracleAttribute()) === missing_oracle
        policy = RSDP.DecimalStringInexact()
        MOI.set(optimizer, RSDP.ExactificationPolicyAttribute(), policy)
        @test MOI.get(optimizer, RSDP.ExactificationPolicyAttribute()) === policy

        for attribute in (
            RSDP.ValidationStatusAttribute(),
            RSDP.CertificateAttribute(),
            RSDP.ValidationReportAttribute(),
            RSDP.OracleResultAttribute(),
            RSDP.DiagnosticsAttribute(),
            RSDP.MaxDenominatorAttribute(),
            RSDP.NumericalOracleAttribute(),
            RSDP.ExactificationPolicyAttribute(),
        )
            @test MOI.supports(optimizer, attribute)
        end

        variable = MOI.add_variable(optimizer)
        @test MOI.is_valid(optimizer, variable)
        MOI.set(optimizer, MOI.VariableName(), variable, "x")
        @test MOI.get(optimizer, MOI.VariableName(), variable) == "x"
        @test MOI.get(optimizer, MOI.VariableIndex, "x") == variable

        equality = MOI.add_constraint(
            optimizer,
            MOI.ScalarAffineFunction(
                [MOI.ScalarAffineTerm(Q(1), variable)],
                Q(0),
            ),
            MOI.EqualTo(Q(1)),
        )
        @test MOI.supports(optimizer, MOI.ConstraintName(), typeof(equality))
        MOI.set(optimizer, MOI.ConstraintName(), equality, "fix_x")
        @test MOI.get(optimizer, MOI.ConstraintName(), equality) == "fix_x"
        @test MOI.get(optimizer, MOI.ConstraintIndex, "fix_x") == equality
        @test MOI.get(optimizer, typeof(equality), "fix_x") == equality
        @test MOI.get(optimizer, MOI.ConstraintSet(), equality) == MOI.EqualTo(Q(1))
        MOI.modify(
            optimizer,
            equality,
            MOI.ScalarCoefficientChange(variable, Q(2)),
        )
        @test only(MOI.get(optimizer, MOI.ConstraintFunction(), equality).terms).coefficient ==
              Q(2)

        @test_throws MOI.UnsupportedConstraint MOI.add_constraint(
            optimizer,
            MOI.ScalarAffineFunction(
                [MOI.ScalarAffineTerm(Q(1), variable)],
                Q(0),
            ),
            MOI.LessThan(Q(1)),
        )

        MOI.delete(optimizer, equality)
        @test !MOI.is_valid(optimizer, equality)
        MOI.empty!(optimizer)
        @test MOI.is_empty(optimizer)

        @test MOI.get(optimizer, RSDP.NumericalOracleAttribute()) isa
              OptimizerMissingOracle
        @test !MOI.supports_constraint(
            optimizer,
            MOI.ScalarAffineFunction{Float64},
            MOI.EqualTo{Float64},
        )
    end

    @testset "extraction and validation failure statuses" begin
        invalid = RSDP.Optimizer()
        variable = MOI.add_variable(invalid)
        MOI.add_constraint(
            invalid.model,
            MOI.ScalarAffineFunction(
                [MOI.ScalarAffineTerm(Q(1), variable)],
                Q(0),
            ),
            MOI.LessThan(Q(1)),
        )
        MOI.optimize!(invalid)
        @test MOI.get(invalid, RSDP.ValidationStatusAttribute()) ==
              RSDP.UNSUPPORTED_MODEL
        @test MOI.get(invalid, MOI.TerminationStatus()) == MOI.INVALID_MODEL
        @test any(contains("unsupported MOI constraint"), RSDP.diagnostics(invalid))

        recovery = RSDP.Optimizer(; oracle = OptimizerRecoveryFailureOracle())
        variable = MOI.add_variable(recovery)
        MOI.add_constraint(
            recovery,
            MOI.ScalarAffineFunction(
                [MOI.ScalarAffineTerm(Q(1), variable)],
                Q(0),
            ),
            MOI.EqualTo(Q(1, 2)),
        )
        MOI.optimize!(recovery)
        @test MOI.get(recovery, MOI.TerminationStatus()) == MOI.OTHER_ERROR
        @test RSDP.validation_status(recovery) ==
              RSDP.RECOVERY_FAILED_DENOMINATOR_LIMIT

        unsupported =
            RSDP.Optimizer(; oracle = OptimizerUnsupportedConeOracle())
        MOI.add_variable(unsupported)
        MOI.optimize!(unsupported)
        @test RSDP.validation_status(unsupported) == RSDP.UNSUPPORTED_CONE
        @test MOI.get(unsupported, MOI.TerminationStatus()) == MOI.INVALID_MODEL

        inexact = RSDP.InexactDataError(
            0.5,
            RSDP.ErrorOnInexact(),
            "test",
            "inexact",
        )
        @test RSDP._extraction_failure_status(inexact) ==
              RSDP.EXACTIFICATION_REQUIRED
        @test RSDP._extraction_failure_status(RSDP.InvalidProblemError("test")) ==
              RSDP.UNSUPPORTED_MODEL
        @test RSDP._extraction_failure_status(ErrorException("test")) ==
              RSDP.UNSUPPORTED_MODEL
    end
end
