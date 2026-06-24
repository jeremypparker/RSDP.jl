using MathOptInterface
using RSDP

const MOI = MathOptInterface
const Q = Rational{BigInt}

model = MOI.Utilities.Model{Q}()
x = MOI.add_variables(model, 2)
MOI.add_constraint(model, MOI.VectorOfVariables(x), MOI.Nonnegatives(2))
MOI.add_constraint(
    model,
    MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(Q[1, 1], x), Q(0)),
    MOI.EqualTo(Q(1)),
)

extracted = extract_moi_problem(model)
@assert extracted.problem.cones !== nothing
