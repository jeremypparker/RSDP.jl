using RSDP
using Test

@testset "RSDP.jl" begin
    include("test_statuses.jl")
    include("test_exactification.jl")
    include("test_exact_affine.jl")
    include("test_cones.jl")
    include("test_psd_checks.jl")
    include("test_adversarial_psd.jl")
    include("test_certificates.jl")
    include("test_rational_recovery.jl")
    include("test_numerical_oracles.jl")
    include("test_hypatia_oracle.jl")
    include("test_docs_examples.jl")
    include("test_spike1.jl")
    include("test_moi_extraction.jl")
    include("test_spike2.jl")
    include("test_boundary.jl")
    include("test_quality.jl")
end
