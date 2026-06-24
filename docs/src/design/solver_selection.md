# Numerical oracle selection

Clarabel is the recommended first hint generator for small standard-cone SDPs because
its interior-point iterates are generally better suited to rational recovery than
first-order solutions. Hypatia is a useful future high-precision and general-cone
adapter. COSMO remains useful for rough large-scale hints.

RSDP accepts an explicit MOI optimizer factory and instantiates it with bridges. In
particular, MOI's PSD dot-scaling bridge keeps irrational off-diagonal scaling outside
the exact problem representation. No solver is loaded or selected implicitly.
