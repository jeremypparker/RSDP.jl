using Test

@testset "exact RREF and nullspace" begin
    Q = RSDP.ExactScalar
    A = Q[0 2 4; 1 1 1; 2 2 2]
    R, pivots = RSDP.exact_rref(A)

    @test pivots == [1, 2]
    @test R == Q[1 0 -1; 0 1 2; 0 0 0]
    N = RSDP.exact_nullspace(A)
    @test N == reshape(Q[1, -2, 1], 3, 1)
    @test A * N == zeros(Q, 3, 1)
    @test RSDP.exact_rank(A) == 2
    @test_throws RSDP.InexactDataError RSDP.exact_rref([1.0 2.0])
    decimal_R, decimal_pivots = RSDP.exact_rref(
        [0.5 1.0];
        policy=RSDP.DecimalStringInexact(),
    )
    @test decimal_R == Q[1 2]
    @test decimal_pivots == [1]
end

@testset "exact affine solve" begin
    Q = RSDP.ExactScalar

    unique_result = RSDP.solve_affine(Q[2 1; 1 -1], Q[5, 1])
    @test RSDP.is_feasible(unique_result)
    @test RSDP.is_unique(unique_result)
    @test unique_result.solution.particular == Q[2, 1]
    @test size(unique_result.solution.nullspace) == (2, 0)

    family_result = RSDP.exact_solve(Q[1 1 1; 2 2 2], Q[3, 6])
    @test RSDP.is_feasible(family_result)
    @test !RSDP.is_unique(family_result)
    @test family_result.solution.pivot_columns == [1]
    @test family_result.solution.free_columns == [2, 3]
    @test family_result.solution.particular == Q[3, 0, 0]
    @test family_result.solution.nullspace == Q[-1 -1; 1 0; 0 1]
    @test RSDP.affine_point(family_result.solution, Q[2, -4]) == Q[5, 2, -4]
    @test RSDP.affine_point(
        family_result.solution,
        [2.0, -4.0];
        policy=RSDP.DecimalStringInexact(),
    ) == Q[5, 2, -4]

    inconsistent_result = RSDP.solve_affine(Q[1 1; 2 2], Q[1, 3])
    @test RSDP.is_infeasible(inconsistent_result)
    @test isnothing(inconsistent_result.solution)
    @test inconsistent_result.diagnostics.rank == 1
    @test inconsistent_result.diagnostics.augmented_rank == 2
    @test_throws RSDP.InconsistentAffineSystemError RSDP.exact_affine_space(
        Q[1 1; 2 2],
        Q[1, 3],
    )
end

@testset "exact conic problem" begin
    Q = RSDP.ExactScalar
    problem = RSDP.ExactConicProblem(Q[1 2; 3 4], Q[5, 11]; c=Q[1, -1])
    @test RSDP.num_constraints(problem) == 2
    @test RSDP.num_variables(problem) == 2
    @test RSDP.solve_affine(problem).solution.particular == Q[1, 2]
    @test RSDP.exact_affine_space(problem).nullspace == zeros(Q, 2, 0)
    @test RSDP.satisfies_affine_constraints(problem, Q[1, 2])

    @test_throws RSDP.InexactDataError RSDP.ExactConicProblem([1.0 0.0], [1.0])
    converted = RSDP.ExactConicProblem(
        [1.0 0.5],
        [2.0];
        policy=RSDP.DecimalStringInexact(),
    )
    @test converted.A == Q[1 1 // 2]
    @test_throws RSDP.InvalidProblemError RSDP.ExactConicProblem(Q[1 2], Q[1, 2])
    @test_throws RSDP.InvalidProblemError RSDP.ExactConicProblem(
        Q[1 2],
        Q[1];
        c=Q[1],
    )
end
