"""
    ExactConicProblem(A, b; c=nothing, policy=ErrorOnInexact())

An exact conic problem represented in flat cone coordinates:

```math
A x = b
```

with optional linear objective `c'x`. Cone-coordinate interpretation is left
to higher-level cone metadata; this core stores `x` as a flat vector.
"""
struct ExactConicProblem
    A::Matrix{ExactScalar}
    b::Vector{ExactScalar}
    c::Union{Nothing,Vector{ExactScalar}}
    cones::Union{Nothing,AbstractConeBlock}
    exactification_policy::AbstractInexactPolicy
    metadata::Dict{Symbol,Any}

    function ExactConicProblem(
        A::Matrix{ExactScalar},
        b::Vector{ExactScalar},
        c::Union{Nothing,Vector{ExactScalar}},
        cones::Union{Nothing,AbstractConeBlock},
        exactification_policy::AbstractInexactPolicy,
        metadata::Dict{Symbol,Any},
    )
        size(A, 1) == length(b) ||
            throw(
                InvalidProblemError(
                    "A has $(size(A, 1)) rows but b has length $(length(b))",
                ),
            )
        isnothing(c) ||
            size(A, 2) == length(c) ||
            throw(
                InvalidProblemError(
                    "A has $(size(A, 2)) columns but c has length $(length(c))",
                ),
            )
        isnothing(cones) ||
            dimension(cones) == size(A, 2) ||
            throw(
                InvalidProblemError(
                    "cone dimension $(dimension(cones)) does not match " *
                    "$(size(A, 2)) variables",
                ),
            )
        new(
            copy(A),
            copy(b),
            isnothing(c) ? nothing : copy(c),
            cones,
            exactification_policy,
            copy(metadata),
        )
    end
end

function ExactConicProblem(
    A::AbstractMatrix,
    b::AbstractVector;
    c=nothing,
    objective=c,
    cones::Union{Nothing,AbstractConeBlock}=nothing,
    policy::AbstractInexactPolicy=DEFAULT_INEXACT_POLICY,
    metadata::AbstractDict{Symbol}=Dict{Symbol,Any}(),
)
    exact_A = exactify(A, policy; context="A")
    exact_b = vec(exactify(b, policy; context="b"))
    exact_c =
        isnothing(objective) ? nothing : vec(exactify(objective, policy; context="c"))
    return ExactConicProblem(
        exact_A,
        exact_b,
        exact_c,
        cones,
        policy,
        Dict{Symbol,Any}(metadata),
    )
end

ExactConicProblem(
    A::AbstractMatrix,
    b::AbstractVector,
    c::AbstractVector;
    policy::AbstractInexactPolicy=DEFAULT_INEXACT_POLICY,
) = ExactConicProblem(A, b; c=c, policy=policy)

function ExactConicProblem(
    A::AbstractMatrix,
    b::AbstractVector,
    cones::AbstractVector{<:AbstractConeBlock};
    objective=nothing,
    c=objective,
    policy::AbstractInexactPolicy=DEFAULT_INEXACT_POLICY,
    metadata::AbstractDict{Symbol}=Dict{Symbol,Any}(),
)
    return ExactConicProblem(
        A,
        b;
        c=c,
        cones=ProductConeBlock(cones),
        policy=policy,
        metadata=metadata,
    )
end

function ExactConicProblem(
    A::AbstractMatrix,
    b::AbstractVector,
    cone::AbstractConeBlock;
    objective=nothing,
    c=objective,
    policy::AbstractInexactPolicy=DEFAULT_INEXACT_POLICY,
    metadata::AbstractDict{Symbol}=Dict{Symbol,Any}(),
)
    return ExactConicProblem(
        A,
        b;
        c=c,
        cones=cone,
        policy=policy,
        metadata=metadata,
    )
end

"""Return the number of scalar equality constraints."""
num_constraints(problem::ExactConicProblem) = size(problem.A, 1)

"""Return the number of flat cone-coordinate variables."""
num_variables(problem::ExactConicProblem) = size(problem.A, 2)

"""Return the exact scalar objective at `x`, or `nothing` for feasibility problems."""
function objective_value(problem::ExactConicProblem, x::AbstractVector)
    isnothing(problem.c) && return nothing
    exact_x = vec(exactify(x; context="x"))
    length(exact_x) == num_variables(problem) ||
        throw(DimensionMismatch("objective point has the wrong dimension"))
    return dot(problem.c, exact_x)
end

"""Return a deterministic SHA-256 fingerprint of exact problem data."""
function problem_hash(problem::ExactConicProblem)
    io = IOBuffer()
    show(io, MIME("text/plain"), problem.A)
    print(io, '\n')
    show(io, MIME("text/plain"), problem.b)
    print(io, '\n')
    show(io, MIME("text/plain"), problem.c)
    print(io, '\n')
    show(io, MIME("text/plain"), problem.cones)
    return bytes2hex(SHA.sha256(take!(io)))
end

"""
    ExactAffineSolution

Exact parameterization `x = particular + nullspace * t` of all solutions to a
consistent affine system. Nullspace basis columns are ordered by increasing
free-variable index.
"""
struct ExactAffineSolution
    particular::Vector{ExactScalar}
    nullspace::Matrix{ExactScalar}
    pivot_columns::Vector{Int}
    free_columns::Vector{Int}

    function ExactAffineSolution(
        particular::Vector{ExactScalar},
        nullspace::Matrix{ExactScalar},
        pivot_columns::AbstractVector{<:Integer},
        free_columns::AbstractVector{<:Integer},
    )
        size(nullspace, 1) == length(particular) ||
            throw(DimensionMismatch("nullspace basis has the wrong row count"))
        size(nullspace, 2) == length(free_columns) ||
            throw(DimensionMismatch("one basis column is required per free variable"))
        new(
            copy(particular),
            copy(nullspace),
            Int[p for p in pivot_columns],
            Int[p for p in free_columns],
        )
    end
end

"""Return the affine dimension of `solution`."""
affine_dimension(solution::ExactAffineSolution) = size(solution.nullspace, 2)

"""Return `true` when an affine solution consists of one point."""
is_unique(solution::ExactAffineSolution) = affine_dimension(solution) == 0

"""
    affine_point(solution[, parameters])

Evaluate the exact affine parameterization. Omitting `parameters` returns the
distinguished particular solution.
"""
function affine_point(
    solution::ExactAffineSolution,
    parameters::AbstractVector=zeros(ExactScalar, affine_dimension(solution)),
    ;
    policy::AbstractInexactPolicy=DEFAULT_INEXACT_POLICY,
)
    length(parameters) == affine_dimension(solution) ||
        throw(
            DimensionMismatch(
                "expected $(affine_dimension(solution)) affine parameters, " *
                "got $(length(parameters))",
            ),
        )
    exact_parameters = exactify(parameters, policy; context="parameters")
    return solution.particular + solution.nullspace * exact_parameters
end

"""
    ExactAffineResult

Status, optional solution parameterization, and structural diagnostics returned
by exact affine solving.
"""
struct ExactAffineResult
    status::SolveStatus
    solution::Union{Nothing,ExactAffineSolution}
    diagnostics::AffineDiagnostics

    function ExactAffineResult(
        status::SolveStatus,
        solution::Union{Nothing,ExactAffineSolution},
        diagnostics::AffineDiagnostics,
    )
        if status == FEASIBLE && isnothing(solution)
            throw(ArgumentError("a feasible result must carry a solution"))
        elseif status == INFEASIBLE && !isnothing(solution)
            throw(ArgumentError("an infeasible result cannot carry a solution"))
        end
        new(status, solution, diagnostics)
    end
end

is_feasible(result::ExactAffineResult) = is_feasible(result.status)
is_infeasible(result::ExactAffineResult) = is_infeasible(result.status)
is_unknown(result::ExactAffineResult) = is_unknown(result.status)

function is_unique(result::ExactAffineResult)
    is_feasible(result) || return false
    return is_unique(something(result.solution))
end
