using Hypatia
using Test
using LinearAlgebra: I

struct UnsupportedOracleCone <: RSDP.AbstractConeBlock end
RSDP.dimension(::UnsupportedOracleCone) = 1

@testset "Hypatia numerical oracle" begin
    Q = RSDP.ExactScalar
    oracle = RSDP.HypatiaOracle()

    nonnegative_problem =
        RSDP.ExactConicProblem(reshape(Q[1], 1, 1), Q[1//2], RSDP.NonnegativeConeBlock(1))
    nonnegative =
        RSDP.validate_with_oracle(nonnegative_problem, oracle; max_denominator = 10^6)
    @test nonnegative.oracle_result.status == RSDP.NUMERICAL_SOLVED_NOT_VALIDATED
    @test nonnegative.certificate !== nothing
    @test nonnegative.certificate.x == Q[1//2]
    @test nonnegative.report.ok

    # The packed order is (1,1), (1,2), (2,2), matching MOI's triangle cone.
    psd_point = Q[1, 1//2, 1]
    psd_problem =
        RSDP.ExactConicProblem(Matrix{Q}(I, 3, 3), psd_point, RSDP.PSDTriangleConeBlock(2))
    psd = RSDP.validate_with_oracle(psd_problem, oracle; max_denominator = 10^6)
    @test psd.certificate !== nothing
    @test psd.certificate.x == psd_point
    @test RSDP.triangle_to_matrix(psd.certificate.x, 2) == Q[1 1//2; 1//2 1]
    @test psd.report.ok

    mixed_point = Q[1//3, 1, 1//3, 1]
    mixed_problem = RSDP.ExactConicProblem(
        Q[
            1 0 0 0
            -1 0 1 0
            0 1 0 0
            0 0 0 1
        ],
        Q[1//3, 0, 1, 1],
        RSDP.AbstractConeBlock[RSDP.NonnegativeConeBlock(1), RSDP.PSDTriangleConeBlock(2)],
    )
    mixed = RSDP.validate_with_oracle(mixed_problem, oracle; max_denominator = 10^6)
    @test mixed.certificate !== nothing
    @test mixed.certificate.x == mixed_point
    @test mixed.report.ok

    unsupported_problem =
        RSDP.ExactConicProblem(reshape(Q[1], 1, 1), Q[0], UnsupportedOracleCone())
    unsupported = RSDP.validate_with_oracle(unsupported_problem, oracle)
    @test isnothing(unsupported.certificate)
    @test unsupported.oracle_result.status == RSDP.UNSUPPORTED_CONE
    @test unsupported.report.status == RSDP.UNSUPPORTED_CONE
    @test any(contains("UnsupportedOracleCone"), unsupported.diagnostics)
end
