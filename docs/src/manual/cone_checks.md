# Cone checks

The zero and nonnegative cones are checked componentwise. PSD triangle vectors use
MathOptInterface's upper-triangle, column-major ordering.

Exact PSD checking uses symmetric elimination with positive diagonal pivots and exact
Schur complements. If every remaining diagonal is zero, the whole remaining matrix
must be zero. This handles singular matrices and zero leading pivots without numerical
eigenvalues.

The initial implementation is intentionally dense and intended for modest blocks.
