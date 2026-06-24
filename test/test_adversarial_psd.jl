module AdversarialPSDTests

using Test

module AdversarialPSDTestImplementation
include(joinpath(@__DIR__, "..", "src", "Cones.jl"))
include(joinpath(@__DIR__, "..", "src", "PSD.jl"))
end

using .AdversarialPSDTestImplementation

adversarial_quadratic_form(matrix, vector) = transpose(vector) * matrix * vector

@testset "adversarial exact PSD checks" begin
    # A rank-one PSD matrix whose only immediately positive diagonal is last.
    permuted_rank_one = [0 0 0; 0 0 0; 0 0 9]
    result = check_psd(permuted_rank_one)
    @test result.is_psd
    @test result.rank == 1
    @test first(result.permutation) == 3

    # The kernel is not aligned with a coordinate axis.
    permuted_kernel = [1 -1 0; -1 1 0; 0 0 2]
    permutation = [3, 1, 2]
    matrix = permuted_kernel[permutation, permutation]
    result = check_psd(matrix)
    @test result.is_psd
    @test result.rank == 2

    # Every original diagonal is positive, but exact elimination exposes
    # indefiniteness. Tiny rational negativity must not be rounded away.
    denominator = big(10)^80
    tiny_negative = Rational{BigInt}[
        1 1
        1 1 - 1//denominator
    ]
    result = check_psd(tiny_negative)
    @test !result.is_psd
    @test result.diagnostic.code == :negative_pivot
    @test adversarial_quadratic_form(tiny_negative, result.witness) < 0

    tiny_positive = Rational{BigInt}[
        1 1
        1 1 + 1//denominator
    ]
    result = check_psd(tiny_positive)
    @test result.is_psd
    @test result.rank == 2

    # A zero diagonal coupled to a positive block becomes a negative Schur
    # pivot after the algorithm safely selects the positive diagonal first.
    delayed_indefinite = Rational{BigInt}[
        0 1 0
        1 2 1
        0 1 2
    ]
    result = check_psd(delayed_indefinite)
    @test !result.is_psd
    @test result.witness !== nothing
    @test adversarial_quadratic_form(delayed_indefinite, result.witness) < 0

    # Large integers exercise BigInt conversion without overflow.
    huge = big(10)^100
    huge_rank_one = [huge^2 huge; huge 1]
    result = check_psd(huge_rank_one)
    @test result.is_psd
    @test result.rank == 1

    # Nonsymmetry is rejected before any elimination can hide it.
    almost_symmetric = Rational{BigInt}[1 1//2; 1//3 1]
    result = check_psd(almost_symmetric)
    @test !result.is_psd
    @test result.diagnostic.code == :nonsymmetric
end

end # module AdversarialPSDTests
