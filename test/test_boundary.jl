using Test

@testset "recovery boundary diagnostics" begin
    Q = RSDP.ExactScalar

    # A tiny numerical cone violation at a uniquely determined boundary point
    # is removed by exact affine reconstruction.
    boundary_problem = RSDP.ExactConicProblem(
        Q[1 0; 0 1],
        Q[0, 1],
        RSDP.NonnegativeConeBlock(2),
    )
    boundary_hint = RSDP.NumericalPrimalHint([-1e-11, 1 + 1e-11])
    recovered = RSDP.recover_primal_certificate(
        boundary_problem,
        boundary_hint;
        max_denominator=10,
        atol=1e-9,
        return_diagnostics=true,
    )
    @test recovered.certificate.x == Q[0, 1]
    @test recovered.diagnostics.max_primal_error <= BigFloat(1e-9)
    @test recovered.diagnostics.affine_dimension == 0
    @test RSDP.check_certificate(boundary_problem, recovered.certificate).ok

    # A genuinely negative exact affine candidate reaches the independent cone
    # boundary and reports that stage, including the checker report.
    outside_problem = RSDP.ExactConicProblem(
        reshape(Q[0], 1, 1),
        Q[0],
        RSDP.NonnegativeConeBlock(1),
    )
    outside_hint = RSDP.NumericalPrimalHint([-1 / 3])
    error = try
        RSDP.recover_primal_certificate(
            outside_problem,
            outside_hint;
            max_denominator=3,
            atol=1e-12,
        )
        nothing
    catch caught
        caught
    end

    @test error isa RSDP.RationalRecoveryError
    @test error.diagnostics.stage == :certificate
    @test error.diagnostics.status == RSDP.RECOVERY_FAILED_CONE
    @test error.diagnostics.certificate_report !== nothing
    @test !error.diagnostics.certificate_report.cones_valid
    @test occursin("cone membership", error.diagnostics.message)
end
