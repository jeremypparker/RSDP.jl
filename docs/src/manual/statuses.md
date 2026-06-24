# Statuses and diagnostics

`ValidationStatus` separates numerical success, exact validation, unsupported input,
recovery failure, and facial-reduction requirements. A `ValidationReport` is the
authoritative result of checking a certificate. Inspect its `diagnostics` rather than
parsing printed text.

In particular, `NUMERICAL_SOLVED_NOT_VALIDATED` is not feasibility,
`VALIDATED_PRIMAL_FEASIBLE` is not optimality, and no infeasibility status is emitted
without an exact dual/Farkas certificate.
