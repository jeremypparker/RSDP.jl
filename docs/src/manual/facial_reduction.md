# Facial reduction

Boundary feasible sets can make rational recovery fragile. RSDP may diagnose that
facial reduction appears necessary, but v0.1 does not automatically change the cone.

A future certified reduction step must store and replay an exact exposing certificate:

```math
s=A^\mathsf{T}y,\qquad b^\mathsf{T}y=0,\qquad s\in K^*,\qquad s\ne0.
```

Heuristic exposing vectors may guide a search but are not part of a validated proof.
For PSD faces, the restricted basis must also be exactly rationally representable.
