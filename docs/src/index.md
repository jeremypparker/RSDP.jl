```@meta
CurrentModule = RSDP
```

# RSDP

RSDP is an exact certificate checker for rational semidefinite programs. A numerical
solver may suggest a point, but only exact affine and cone checks validate it.

The current package proves primal feasibility and evaluates the attained rational
objective. It does not turn primal feasibility into an optimality claim.

## Core guarantees

- Model coefficients are exact or explicitly exactified.
- Certificate checking uses exact arithmetic and no numerical solver.
- Singular PSD matrices are handled without floating-point eigenvalues.
- Certificate problem hashes detect accidental or malicious reuse.
- Unsupported claims fail with structured statuses and diagnostics.

Start with [Quick start](@ref) and read [Limitations](@ref) before relying on a
certificate in a proof pipeline.

The tutorials show a [JuMP rational SDP](@ref) and an
[ordinary sum-of-squares example](@ref "Ordinary sum of squares") using
Hypatia only as a numerical oracle.
