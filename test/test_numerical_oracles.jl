using Test

struct SuccessfulFakeOracle <: RSDP.AbstractNumericalOracle end
struct MissingPrimalOracle <: RSDP.AbstractNumericalOracle end
struct WrongDimensionOracle <: RSDP.AbstractNumericalOracle end
struct RecoveryFailureOracle <: RSDP.AbstractNumericalOracle end

function RSDP.solve_oracle(::RSDP.ExactConicProblem, ::SuccessfulFakeOracle; kwargs...)
    return RSDP.NumericalOracleResult(
        RSDP.NUMERICAL_SOLVED_NOT_VALIDATED;
        primal = [0.5],
        diagnostics = ["fake numerical hint"],
    )
end

function RSDP.solve_oracle(::RSDP.ExactConicProblem, ::MissingPrimalOracle; kwargs...)
    return RSDP.NumericalOracleResult(
        RSDP.NUMERICAL_ORACLE_FAILED;
        diagnostics = ["fake solver produced no point"],
    )
end

function RSDP.solve_oracle(::RSDP.ExactConicProblem, ::WrongDimensionOracle; kwargs...)
    return RSDP.NumericalOracleResult(
        RSDP.NUMERICAL_SOLVED_NOT_VALIDATED;
        primal = [0.5, 0.5],
    )
end

function RSDP.solve_oracle(::RSDP.ExactConicProblem, ::RecoveryFailureOracle; kwargs...)
    return RSDP.NumericalOracleResult(RSDP.NUMERICAL_SOLVED_NOT_VALIDATED; primal = [-0.5])
end

@testset "numerical oracle validation orchestration" begin
    Q = RSDP.ExactScalar
    problem =
        RSDP.ExactConicProblem(reshape(Q[1], 1, 1), Q[1//2], RSDP.NonnegativeConeBlock(1))

    success = RSDP.validate_with_oracle(problem, SuccessfulFakeOracle())
    @test success.certificate !== nothing
    @test success.certificate.x == Q[1//2]
    @test success.report.ok
    @test success.report.status == RSDP.VALIDATED_PRIMAL_FEASIBLE

    missing = RSDP.validate_with_oracle(problem, MissingPrimalOracle())
    @test isnothing(missing.certificate)
    @test !missing.report.ok
    @test missing.report.status == RSDP.NUMERICAL_ORACLE_FAILED
    @test any(contains("no point"), missing.diagnostics)

    wrong = RSDP.validate_with_oracle(problem, WrongDimensionOracle())
    @test isnothing(wrong.certificate)
    @test !wrong.report.ok
    @test any(contains("expected 1"), wrong.diagnostics)

    recovery =
        RSDP.validate_with_oracle(problem, RecoveryFailureOracle(); atol = 0, rtol = 0)
    @test isnothing(recovery.certificate)
    @test !recovery.report.ok
    @test recovery.report.status == RSDP.RECOVERY_FAILED_DENOMINATOR_LIMIT
    @test any(contains("rational recovery failed"), recovery.diagnostics)
end
