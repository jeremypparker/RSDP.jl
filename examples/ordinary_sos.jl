# Ordinary SOS validation in v0.1 uses the bridged MOI SDP. Build the
# SumOfSquares model in an integration environment, bridge it to affine
# equalities and PSD triangle cones, then call `RSDP.extract_moi_problem`.
#
# Human-readable polynomial identity reconstruction is intentionally separate
# from conic certificate validation.
using RSDP

println("See docs/src/design/weighted_sos.md for the supported bridged workflow.")
