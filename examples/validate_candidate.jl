using RSDP

const Q = Rational{BigInt}
problem = ExactConicProblem(Q[1 1], Q[1], NonnegativeConeBlock(2); objective = Q[1, 2])

certificate = make_primal_certificate(problem, Q[1//3, 2//3])
report = check_certificate(problem, certificate)

@assert report.status == VALIDATED_PRIMAL_FEASIBLE
@assert report.computed_objective == 5 // 3
