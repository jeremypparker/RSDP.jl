using Test

@testset "cone blocks and exact membership" begin
    zero = ZeroConeBlock(3)
    nonnegative = NonnegativeConeBlock(2)
    psd = PSDTriangleConeBlock(3)
    product = ProductConeBlock(zero, nonnegative, psd)

    @test dimension(zero) == 3
    @test dimension(nonnegative) == 2
    @test dimension(psd) == 6
    @test matrix_dimension(psd) == 3
    @test ambient_dimension(product) == 11
    @test length(ProductConeBlock()) == 0
    @test_throws ArgumentError ZeroConeBlock(-1)
    @test_throws ArgumentError PSDTriangleConeBlock(-1)

    @test triangle_dimension(0) == 0
    @test triangle_dimension(4) == 10
    @test triangle_side_dimension(10) == 4
    @test_throws ArgumentError triangle_side_dimension(5)

    packed = [1, 2 // 3, 4, -5, 6, 7]
    matrix = triangle_to_matrix(packed, 3)
    @test eltype(matrix) == Rational{BigInt}
    @test matrix == Rational{BigInt}[
        1 2//3 -5
        2//3 4 6
        -5 6 7
    ]
    @test matrix_to_triangle(matrix) == Rational{BigInt}.(packed)
    @test reconstruct_triangle(packed) == matrix
    @test_throws DimensionMismatch triangle_to_matrix([1, 2], 2)
    @test_throws ArgumentError triangle_to_matrix([1.0], 1)
    @test_throws ArgumentError matrix_to_triangle([1 2; 3 4])

    @test check_cone_membership(zero, [0, 0 // 2, 0]).is_member
    zero_failure = check_cone_membership(zero, [0, 1, 0])
    @test !zero_failure.is_member
    @test only(zero_failure.diagnostics).code == :nonzero_in_zero_cone

    @test in_cone(nonnegative, [0, 3 // 2])
    nonnegative_failure = check_cone(nonnegative, [0, -1])
    @test !nonnegative_failure.is_member
    @test only(nonnegative_failure.diagnostics).index == 2

    psd_two = PSDTriangleConeBlock(2)
    @test in_cone(psd_two, [1, 1, 1])
    psd_failure = check_cone_membership(psd_two, [1, 2, 1])
    @test !psd_failure.is_member
    @test only(psd_failure.diagnostics).code == :not_positive_semidefinite
    @test only(psd_failure.diagnostics).psd_result.witness !== nothing

    nested = ProductConeBlock(ZeroConeBlock(1), ProductConeBlock(nonnegative, psd_two))
    nested_failure = check_cone_membership(nested, [0, 1, -1, 1, 2, 1])
    @test !nested_failure.is_member
    @test [diagnostic.path for diagnostic in nested_failure.diagnostics] == [[2, 1], [2, 2]]
    @test check_cones((ZeroConeBlock(1), NonnegativeConeBlock(1)), [0, 2]).is_member

    mismatch = check_cone_membership(product, zeros(Int, 2))
    @test !mismatch.is_member
    @test only(mismatch.diagnostics).code == :dimension_mismatch

    inexact = check_cone_membership(NonnegativeConeBlock(1), [0.0])
    @test !inexact.is_member
    @test only(inexact.diagnostics).code == :inexact_input
end
