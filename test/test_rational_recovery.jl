using Test

@testset "affine-coordinate rational recovery" begin
    Q = RSDP.ExactScalar
    problem = RSDP.ExactConicProblem(
        Q[1 1],
        Q[1],
        RSDP.NonnegativeConeBlock(2),
        objective=Q[1, 2],
    )
    hint = RSDP.NumericalPrimalHint(
        [0.333333333333, 0.666666666667];
        objective=1.666666666667,
    )

    options = RSDP.RecoveryOptions(
        max_denominator=100,
        atol=1e-9,
        rtol=1e-9,
    )
    recovered = RSDP.recover_primal_certificate(
        problem,
        hint,
        options;
        return_diagnostics=true,
    )

    @test recovered.certificate.x == Q[1 // 3, 2 // 3]
    @test recovered.certificate.objective == 5 // 3
    @test recovered.diagnostics.stage == :success
    @test recovered.diagnostics.status == RSDP.VALIDATED_PRIMAL_FEASIBLE
    @test recovered.diagnostics.affine_dimension == 1
    @test recovered.diagnostics.certificate_report.valid
    @test RSDP.check_certificate(problem, recovered.certificate).ok
    @test RSDP.recover_primal_certificate(problem, hint; options=options).x ==
          Q[1 // 3, 2 // 3]

    poor_bound = RSDP.NumericalPrimalHint([1 / 3, 2 / 3])
    error = try
        RSDP.recover_primal_certificate(
            problem,
            poor_bound;
            max_denominator=2,
            atol=1e-5,
            rtol=0,
        )
        nothing
    catch caught
        caught
    end
    @test error isa RSDP.RationalRecoveryError
    @test error.diagnostics.stage == :proximity
    @test error.diagnostics.status == RSDP.RECOVERY_FAILED_DENOMINATOR_LIMIT
    @test occursin("bounded-denominator", error.diagnostics.message)

    wrong_objective =
        RSDP.NumericalPrimalHint([1 / 3, 2 / 3]; objective=9.0)
    objective_error = try
        RSDP.recover_primal_certificate(
            problem,
            wrong_objective;
            max_denominator=100,
            atol=1e-8,
        )
        nothing
    catch caught
        caught
    end
    @test objective_error isa RSDP.RationalRecoveryError
    @test objective_error.diagnostics.stage == :objective
    @test objective_error.diagnostics.status == RSDP.RECOVERY_FAILED
end
