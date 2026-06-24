using Test
using LinearAlgebra: dot
import MathOptInterface as MOI

@testset "MOI extraction spike: fixed 2x2 PSD matrix" begin
    Q = Rational{BigInt}
    model = MOI.Utilities.Model{Q}()
    diagonal_1, off_diagonal, diagonal_2 = MOI.add_variables(model, 3)

    MOI.add_constraint(
        model,
        MOI.VectorAffineFunction{Q}(
            [
                MOI.VectorAffineTerm(
                    1,
                    MOI.ScalarAffineTerm(Q(1), diagonal_1),
                ),
                MOI.VectorAffineTerm(
                    2,
                    MOI.ScalarAffineTerm(Q(1), diagonal_2),
                ),
            ],
            Q[-1, -1],
        ),
        MOI.Zeros(2),
    )
    MOI.add_constraint(
        model,
        MOI.VectorOfVariables([diagonal_1, off_diagonal, diagonal_2]),
        MOI.PositiveSemidefiniteConeTriangle(2),
    )
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Q}}(),
        MOI.ScalarAffineFunction(
            [MOI.ScalarAffineTerm(Q(1), off_diagonal)],
            Q(0),
        ),
    )

    extracted = RSDP.extract_moi(model)
    @test extracted.block_columns == [1:3, 4:6, 7:9]
    @test extracted.blocks == [
        RSDP.NonnegativeConeBlock(3),
        RSDP.NonnegativeConeBlock(3),
        RSDP.PSDTriangleConeBlock(2),
    ]

    # u = (1, 0, 1), represented by positive-minus-negative parts, and
    # y = (1, 0, 1) is the corresponding packed PSD matrix.
    coordinates = Q[1, 0, 1, 0, 0, 0, 1, 0, 1]
    @test extracted.A * coordinates == extracted.b
    @test RSDP.in_cone(extracted.cone, coordinates)
    @test RSDP.recover_moi_variables(extracted, coordinates) == Q[1, 0, 1]
    @test dot(extracted.c, coordinates) + extracted.objective_constant == Q(0)
end
