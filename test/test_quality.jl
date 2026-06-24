using Aqua
using JET

@testset "quality checks" begin
    Aqua.test_all(RSDP; ambiguities=false)
    JET.test_package(RSDP; target_modules=(RSDP,))
end
