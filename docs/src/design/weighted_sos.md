# Weighted SOS

The intended human-readable certificate is

```math
p=\sum_i w_i\,b_i^\mathsf{T}Q_i b_i,\qquad Q_i\succeq0.
```

There are two validation levels:

1. Validate the exact bridged MOI SDP. This proves the conic Gram certificate.
2. When stable metadata are available, also verify exact weights, bases, Gram blocks,
   and the polynomial identity.

RSDP does not infer lost polynomial metadata from an arbitrary bridged model. Any
SumOfSquares-specific implementation is isolated from the core checker.

## SumOfSquares.jl v0.8.0 findings

`WeightedSOSCone` publicly exposes `basis`, `gram_bases`, and `weights`, which is
enough to define the exact symbolic map. Numerical Gram matrices are reliably
available through `GramMatrixAttribute(multiplier_index=i)` on the KernelBridge
route, using one-based direct block indices.

In v0.8.0, ImageBridge has known limitations for non-unit weights and multiple Gram
blocks, while LowRankBridge does not expose `GramMatrixAttribute`. RSDP therefore
documents Kernel-compatible bridged validation as the supported early path and
returns an unsupported result for routes that cannot provide a complete witness.
