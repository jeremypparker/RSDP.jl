# SumOfSquares.jl v0.8 WeightedSOSCone support is version-sensitive.
# The KernelBridge route exposes one GramMatrixAttribute per gram basis and is
# suitable for numerical hints. ImageBridge and LowRankBridge do not expose a
# uniform public Gram-certificate path in v0.8.0.
using RSDP

println("Validate the bridged MOI model; see docs/src/design/weighted_sos.md.")
