# Scaling

Dense rational nullspaces and Schur complements can grow quickly in both memory and
coefficient size. RSDP records dimensions and emits diagnostics before known
high-growth operations.

The design preserves cone blocks so future implementations can add block-decomposed
recovery, sparse fraction-free elimination, modular rank and nullspace computation,
and rational reconstruction without changing certificate semantics.

Supplied-candidate validation avoids nullspace construction and is the preferred path
for large problems when an exact point is already available.
