using Test

@testset "exactification policies" begin
    Q = RSDP.ExactScalar

    @test RSDP.exactify(3) == Q(3)
    @test RSDP.exactify(2 // 6) == Q(1, 3)
    @test RSDP.exactify("-.125e1") == Q(-5, 4)
    @test RSDP.exactify("12/18") == Q(2, 3)
    @test_throws RSDP.InexactDataError RSDP.exactify(0.5)
    @test_throws RSDP.InexactDataError RSDP.exactify(Inf, RSDP.RationalizeInexact())
    @test_throws RSDP.InexactDataError RSDP.exactify(NaN, RSDP.DecimalStringInexact())

    @test RSDP.exactify(0.1, RSDP.RationalizeInexact()) == Q(1, 10)
    @test RSDP.exactify(0.1, RSDP.RationalizeInexact(tol=0.0)) ==
          Q(BigInt(3602879701896397), BigInt(36028797018963968))
    @test RSDP.exactify(0.1, RSDP.RationalizeInexact(0.0)) ==
          Q(BigInt(3602879701896397), BigInt(36028797018963968))
    @test RSDP.exactify(0.1, RSDP.DecimalStringInexact()) == Q(1, 10)
    @test RSDP.exactify(nextfloat(0.1), RSDP.DecimalStringInexact()) ==
          Q(5000000000000001, 50000000000000000)

    exact_matrix = RSDP.exactify(
        [0.5 1.25; -2.0 3.0],
        RSDP.DecimalStringInexact();
        context="A",
    )
    @test exact_matrix == Q[1 // 2 5 // 4; -2 3]
    @test eltype(exact_matrix) == Q
    @test_throws ArgumentError RSDP.RationalizeInexact(tolerance=1e-6, tol=1e-6)
end
