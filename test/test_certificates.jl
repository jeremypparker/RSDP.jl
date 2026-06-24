using Test

@testset "exact primal certificates" begin
    Q = RSDP.ExactScalar
    problem = RSDP.ExactConicProblem(
        Q[1 1],
        Q[1],
        RSDP.NonnegativeConeBlock(2),
        objective = Q[1, 2],
        metadata = Dict(:name => "certificate test"),
    )
    x = Q[1//3, 2//3]

    certificate = RSDP.make_primal_certificate(problem, x)
    @test certificate isa RSDP.ExactPrimalCertificate
    @test certificate.certificate_version == RSDP.CERTIFICATE_VERSION
    @test certificate.exactification_policy isa RSDP.ErrorOnInexact
    @test certificate.metadata[:name] == "certificate test"
    @test eltype(certificate.x) == Rational{BigInt}
    @test certificate.x_exact == x
    @test certificate.objective == 5 // 3
    @test RSDP.check_certificate(problem, certificate).ok
    @test RSDP.check_certificate(problem, certificate; diagnostics = false)

    report = RSDP.check_certificate(problem, certificate; diagnostics = true)
    @test report isa RSDP.ValidationReport
    @test report.valid
    @test report.status == RSDP.VALIDATED_PRIMAL_FEASIBLE
    @test report.version_valid
    @test report.problem_hash_valid
    @test report.affine_valid
    @test report.cones_valid
    @test report.objective_valid
    @test isempty(report.diagnostics)

    modified_problem = RSDP.ExactConicProblem(
        Q[1 1],
        Q[2],
        RSDP.NonnegativeConeBlock(2),
        objective = Q[1, 2],
    )
    modified_report =
        RSDP.check_certificate(modified_problem, certificate; diagnostics = true)
    @test !modified_report.valid
    @test !modified_report.problem_hash_valid
    @test !modified_report.affine_valid

    affine_failure =
        RSDP.PrimalCertificate(certificate.problem_hash, Q[1//2, 1//3], Q(7 // 6))
    affine_report = RSDP.check_certificate(problem, affine_failure; diagnostics = true)
    @test !affine_report.affine_valid

    cone_failure = RSDP.make_primal_certificate(problem, Q[2, -1])
    cone_report = RSDP.check_certificate(problem, cone_failure; diagnostics = true)
    @test cone_report.affine_valid
    @test !cone_report.cones_valid
    @test any(contains("cone membership"), cone_report.diagnostics)

    objective_failure = RSDP.PrimalCertificate(certificate.problem_hash, certificate.x, 0)
    objective_report =
        RSDP.check_certificate(problem, objective_failure; diagnostics = true)
    @test !objective_report.objective_valid
    @test objective_report.computed_objective == 5 // 3

    wrong_version = RSDP.PrimalCertificate(
        certificate.problem_hash,
        certificate.x,
        certificate.objective;
        certificate_version = v"99.0.0",
    )
    version_report = RSDP.check_certificate(problem, wrong_version)
    @test !version_report.ok
    @test !version_report.version_valid
    @test version_report.status == RSDP.CERTIFICATE_CHECK_FAILED

    no_objective = RSDP.make_primal_certificate(problem, x; include_objective = false)
    @test isnothing(no_objective.objective)
    @test RSDP.check_certificate(problem, no_objective).ok

    affine_only = RSDP.ExactConicProblem(reshape(Q[1], 1, 1), Q[1])
    @test RSDP.check_certificate(
        affine_only,
        RSDP.make_primal_certificate(affine_only, Q[1]),
    ).ok
end
