# Contributing to RSDP

RSDP treats numerical results as hints and exact certificate checks as proofs. Changes to
recovery code must not weaken the independent checker.

Every user-facing change should include deterministic tests, adversarial cases where
appropriate, docstrings, and documentation. New floating-point conversion paths must be
explicitly selected by the caller. New infeasibility or optimality statuses require exact
dual certificate checks before they can be emitted.

Optional solver and polynomial-system integrations belong behind package extensions or in
integration-test environments. The core package should remain usable without a numerical
solver.
