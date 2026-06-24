using Aqua
using JET

@testset "quality checks" begin
    Aqua.test_all(RSDP; ambiguities = false)
    if VERSION >= v"1.12.0" && VERSION < v"1.13.0-"
        using JET
        JET.test_package(RSDP; target_modules = (RSDP,))
    end
end
