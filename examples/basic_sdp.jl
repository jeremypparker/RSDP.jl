using RSDP

const Q = Rational{BigInt}

# x[1:2] is nonnegative and x[3:5] packs [1 0; 0 1].
problem = ExactConicProblem(
    Q[1 1 0 0 0],
    Q[1],
    AbstractConeBlock[
        NonnegativeConeBlock(2),
        PSDTriangleConeBlock(2),
    ],
)

point = Q[1 // 2, 1 // 2, 1, 0, 1]
certificate = make_primal_certificate(problem, point)
report = check_certificate(problem, certificate)
@assert report.ok
