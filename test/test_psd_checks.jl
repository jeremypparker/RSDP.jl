module PSDChecksTests

using Test

module PSDTestImplementation
include(joinpath(@__DIR__, "..", "src", "Cones.jl"))
include(joinpath(@__DIR__, "..", "src", "PSD.jl"))
end

using .PSDTestImplementation

quadratic_form(matrix, vector) = transpose(vector) * matrix * vector

@testset "exact PSD checks" begin
    positive_definite = [2 1; 1 2]
    result = check_psd(positive_definite)
    @test result.is_psd
    @test result.rank == 2
    @test result.diagnostic.code == :positive_definite
    @test result.pivots == Rational{BigInt}[2, 3 // 2]
    @test result.witness === nothing

    singular = [1 1 0; 1 1 0; 0 0 0]
    result = exact_psd_check(singular)
    @test result.is_psd
    @test result.rank == 1
    @test result.diagnostic.code == :positive_semidefinite
    @test count(iszero, result.pivots) == 2

    zero_leading_pivot = [0 0; 0 3]
    result = check_psd(zero_leading_pivot)
    @test result.is_psd
    @test result.rank == 1
    @test result.permutation == [2, 1]

    indefinite = Rational{BigInt}[1 2; 2 1]
    result = check_psd(indefinite)
    @test !result.is_psd
    @test result.diagnostic.code == :negative_pivot
    @test result.witness !== nothing
    @test quadratic_form(indefinite, result.witness) < 0
    @test !is_psd_exact(indefinite)
    @test !is_psd(indefinite)

    zero_diagonal = Rational{BigInt}[0 1; 1 0]
    result = check_psd(zero_diagonal)
    @test !result.is_psd
    @test result.diagnostic.code == :zero_diagonal_nonzero_row
    @test quadratic_form(zero_diagonal, result.witness) < 0

    nonsymmetric = check_psd([1 2; 3 4])
    @test !nonsymmetric.is_psd
    @test nonsymmetric.diagnostic.code == :nonsymmetric

    nonsquare = check_psd(ones(Int, 2, 3))
    @test !nonsquare.is_psd
    @test nonsquare.diagnostic.code == :nonsquare

    inexact = check_psd([1.0 0.0; 0.0 1.0])
    @test !inexact.is_psd
    @test inexact.diagnostic.code == :inexact_input

    empty_result = check_psd(zeros(Int, 0, 0))
    @test empty_result.is_psd
    @test empty_result.rank == 0
    @test empty_result.diagnostic.code == :positive_definite
end

end # module PSDChecksTests
