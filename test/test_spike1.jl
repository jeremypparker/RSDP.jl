using Test

@testset "spike-scale affine recovery" begin
    Q = RSDP.ExactScalar
    problem = RSDP.ExactConicProblem(
        Q[-1000 1 0; 1 0 1],
        Q[0, 1],
        RSDP.NonnegativeConeBlock(3),
        objective = Q[0, 0, 1],
    )

    # The three entries carry independent noise at very different scales.
    # Rationalizing them independently need not preserve either exact equation.
    hint = RSDP.NumericalPrimalHint(
        [Float64(1 // 7) + 1e-7, Float64(1000 // 7) + 2e-6, Float64(6 // 7) - 1e-7];
        objective = Float64(6 // 7),
    )

    certificate = RSDP.recover_primal_certificate(
        problem,
        hint;
        max_denominator = 20,
        atol = 3e-6,
        rtol = 0,
    )
    @test certificate.x == Q[1//7, 1000//7, 6//7]
    @test problem.A * certificate.x == problem.b
    @test RSDP.check_certificate(problem, certificate).ok
end
