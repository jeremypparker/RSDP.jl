using BenchmarkTools
using MathOptInterface
using RSDP

const MOI = MathOptInterface
const MODEL = let
    model = MOI.Utilities.Model{Rational{BigInt}}()
    x = MOI.add_variables(model, 20)
    MOI.add_constraint(model, MOI.VectorOfVariables(x), MOI.Nonnegatives(20))
    model
end

@benchmark extract_moi_problem($MODEL)
