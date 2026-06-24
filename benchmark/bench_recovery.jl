using BenchmarkTools
using RSDP

const Q = Rational{BigInt}
const RECOVERY_PROBLEM = ExactConicProblem(
    Q[1 1 0; 0 0 1],
    Q[1, 1],
    AbstractConeBlock[NonnegativeConeBlock(2), NonnegativeConeBlock(1)],
)
const RECOVERY_HINT = NumericalPrimalHint([0.5000000001, 0.4999999999, 1.0])

@benchmark recover_primal_certificate($RECOVERY_PROBLEM, $RECOVERY_HINT)
