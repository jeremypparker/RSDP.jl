using Test

@testset "documentation examples" begin
    jump_result = include(joinpath(@__DIR__, "..", "examples", "jump_rational_sdp.jl"))
    @test jump_result.certificate !== nothing
    @test jump_result.report.ok

    sos_result = include(joinpath(@__DIR__, "..", "examples", "ordinary_sos.jl"))
    @test sos_result.certificate !== nothing
    @test sos_result.report.ok
end
