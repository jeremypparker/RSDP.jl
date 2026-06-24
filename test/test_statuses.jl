using Test

@testset "exact affine statuses" begin
    @test RSDP.is_unknown(RSDP.UNKNOWN)
    @test RSDP.is_feasible(RSDP.FEASIBLE)
    @test RSDP.is_infeasible(RSDP.INFEASIBLE)
    @test !RSDP.is_feasible(RSDP.UNKNOWN)
    @test RSDP.ExactAffineStatus === RSDP.SolveStatus
end
