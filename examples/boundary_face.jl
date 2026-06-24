using RSDP
using LinearAlgebra

const Q = Rational{BigInt}

# The exact affine equations force the PSD matrix diag(0, 1), which lies on a
# proper face. Recovery succeeds because the face is already encoded exactly.
problem = ExactConicProblem(Matrix{Q}(I, 3, 3), Q[0, 0, 1], PSDTriangleConeBlock(2))
hint = NumericalPrimalHint([-1e-12, 2e-12, 1 + 1e-12])
certificate = recover_primal_certificate(problem, hint; atol = 1e-9)
@assert check_certificate(problem, certificate).ok
