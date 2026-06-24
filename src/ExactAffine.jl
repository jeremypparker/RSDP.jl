"""
    solve_affine(A, b; policy=ErrorOnInexact())
    solve_affine(problem)

Solve `A*x = b` exactly and return an [`ExactAffineResult`](@ref). A feasible
result parameterizes every solution as a particular point plus a deterministic
nullspace basis.
"""
function solve_affine(
    A::AbstractMatrix,
    b::AbstractVector;
    policy::AbstractInexactPolicy = DEFAULT_INEXACT_POLICY,
)
    data = exact_linear_solve(A, b; policy = policy)
    diagnostics = AffineDiagnostics(
        size(A, 1),
        size(A, 2),
        data.rank,
        data.augmented_rank,
        data.pivot_columns,
        data.free_columns,
    )

    if !data.consistent
        return ExactAffineResult(INFEASIBLE, nothing, diagnostics)
    end

    solution = ExactAffineSolution(
        something(data.particular),
        data.nullspace,
        data.pivot_columns,
        data.free_columns,
    )
    return ExactAffineResult(FEASIBLE, solution, diagnostics)
end

solve_affine(A::AbstractMatrix, b::AbstractVector, policy::AbstractInexactPolicy) =
    solve_affine(A, b; policy = policy)

solve_affine(problem::ExactConicProblem) = solve_affine(problem.A, problem.b)

"""
Alias for [`solve_affine`](@ref).
"""
exact_solve(args...; kwargs...) = solve_affine(args...; kwargs...)

"""
Alias for [`solve_affine`](@ref).
"""
solve_exact(args...; kwargs...) = solve_affine(args...; kwargs...)

"""
    exact_affine_space(A, b; policy=ErrorOnInexact())
    exact_affine_space(problem)

Return the exact affine-space parameterization for a consistent system.
Unlike [`solve_affine`](@ref), this convenience API throws
[`InconsistentAffineSystemError`](@ref) when no affine space exists.
"""
function exact_affine_space(
    A::AbstractMatrix,
    b::AbstractVector;
    policy::AbstractInexactPolicy = DEFAULT_INEXACT_POLICY,
)
    result = solve_affine(A, b; policy = policy)
    if is_infeasible(result)
        throw(
            InconsistentAffineSystemError(
                result.diagnostics.rank,
                result.diagnostics.augmented_rank,
            ),
        )
    end
    return something(result.solution)
end

exact_affine_space(A::AbstractMatrix, b::AbstractVector, policy::AbstractInexactPolicy) =
    exact_affine_space(A, b; policy = policy)

exact_affine_space(problem::ExactConicProblem) = exact_affine_space(problem.A, problem.b)

"""
    satisfies_affine_constraints(A, b, x; policy=ErrorOnInexact())

Return whether `x` satisfies `A*x == b` after exactification under `policy`.
"""
function satisfies_affine_constraints(
    A::AbstractMatrix,
    b::AbstractVector,
    x::AbstractVector;
    policy::AbstractInexactPolicy = DEFAULT_INEXACT_POLICY,
)
    size(A, 2) == length(x) || throw(
        DimensionMismatch("A has $(size(A, 2)) columns but x has length $(length(x))"),
    )
    size(A, 1) == length(b) ||
        throw(DimensionMismatch("A has $(size(A, 1)) rows but b has length $(length(b))"))
    exact_A = exactify(A, policy; context = "A")
    exact_b = vec(exactify(b, policy; context = "b"))
    exact_x = vec(exactify(x, policy; context = "x"))
    return exact_A * exact_x == exact_b
end

satisfies_affine_constraints(
    problem::ExactConicProblem,
    x::AbstractVector;
    policy::AbstractInexactPolicy = DEFAULT_INEXACT_POLICY,
) = satisfies_affine_constraints(problem.A, problem.b, x; policy = policy)

"""
Alias for [`satisfies_affine_constraints`](@ref).
"""
check_affine(args...; kwargs...) = satisfies_affine_constraints(args...; kwargs...)
