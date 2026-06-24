"""
    exact_rref(A; policy=ErrorOnInexact())

Compute the reduced row-echelon form of `A` over [`ExactScalar`](@ref).
Returns `(R, pivot_columns)`. Pivot selection is deterministic: columns are
visited left-to-right and the topmost available nonzero pivot is chosen.
"""
function exact_rref(
    A::AbstractMatrix;
    policy::AbstractInexactPolicy=DEFAULT_INEXACT_POLICY,
)
    R = exactify(A, policy; context="A")
    pivots = _exact_rref!(R)
    return R, pivots
end

exact_rref(A::AbstractMatrix, policy::AbstractInexactPolicy) =
    exact_rref(A; policy=policy)

"""Alias for [`exact_rref`](@ref)."""
rref_exact(args...; kwargs...) = exact_rref(args...; kwargs...)

function _exact_rref!(R::Matrix{ExactScalar})
    row_count, column_count = size(R)
    pivots = Int[]
    pivot_row = 1

    for column in 1:column_count
        pivot_row > row_count && break
        selected_row = findfirst(row -> !iszero(R[row, column]), pivot_row:row_count)
        isnothing(selected_row) && continue
        selected_row = pivot_row + selected_row - 1

        if selected_row != pivot_row
            for j in 1:column_count
                R[pivot_row, j], R[selected_row, j] =
                    R[selected_row, j], R[pivot_row, j]
            end
        end

        pivot_value = R[pivot_row, column]
        for j in 1:column_count
            R[pivot_row, j] /= pivot_value
        end

        for row in 1:row_count
            row == pivot_row && continue
            factor = R[row, column]
            iszero(factor) && continue
            for j in 1:column_count
                R[row, j] -= factor * R[pivot_row, j]
            end
        end

        push!(pivots, column)
        pivot_row += 1
    end
    return pivots
end

"""Return the exact rank of `A`."""
exact_rank(A::AbstractMatrix; policy::AbstractInexactPolicy=DEFAULT_INEXACT_POLICY) =
    length(last(exact_rref(A; policy=policy)))

"""
    exact_nullspace(A; policy=ErrorOnInexact())

Return a deterministic exact basis for `null(A)` as matrix columns. Free
variables are visited in increasing column order and each corresponding basis
vector sets that free variable to one.
"""
function exact_nullspace(
    A::AbstractMatrix;
    policy::AbstractInexactPolicy=DEFAULT_INEXACT_POLICY,
)
    R, pivots = exact_rref(A; policy=policy)
    return _nullspace_from_rref(R, pivots, size(A, 2))
end

exact_nullspace(A::AbstractMatrix, policy::AbstractInexactPolicy) =
    exact_nullspace(A; policy=policy)

"""Alias for [`exact_nullspace`](@ref)."""
nullspace_exact(args...; kwargs...) = exact_nullspace(args...; kwargs...)

function _nullspace_from_rref(
    R::Matrix{ExactScalar},
    pivots::Vector{Int},
    variable_count::Int,
)
    coefficient_pivots = Int[p for p in pivots if p <= variable_count]
    pivot_set = Set(coefficient_pivots)
    free_columns = Int[j for j in 1:variable_count if !(j in pivot_set)]
    basis = zeros(ExactScalar, variable_count, length(free_columns))

    for (basis_column, free_column) in enumerate(free_columns)
        basis[free_column, basis_column] = one(ExactScalar)
        for (pivot_row, pivot_column) in enumerate(coefficient_pivots)
            basis[pivot_column, basis_column] = -R[pivot_row, free_column]
        end
    end
    return basis
end

"""
    exact_linear_solve(A, b; policy=ErrorOnInexact())

Perform exact Gauss-Jordan elimination on `[A b]`. The returned named tuple
contains consistency, a deterministic particular solution (or `nothing`),
nullspace basis, ranks, pivot/free columns, and the augmented RREF.
"""
function exact_linear_solve(
    A::AbstractMatrix,
    b::AbstractVector;
    policy::AbstractInexactPolicy=DEFAULT_INEXACT_POLICY,
)
    size(A, 1) == length(b) ||
        throw(
            DimensionMismatch(
                "A has $(size(A, 1)) rows but b has length $(length(b))",
            ),
        )
    exact_A = exactify(A, policy; context="A")
    exact_b = vec(exactify(b, policy; context="b"))
    variable_count = size(exact_A, 2)
    augmented = hcat(exact_A, exact_b)
    all_pivots = _exact_rref!(augmented)
    coefficient_pivots = Int[p for p in all_pivots if p <= variable_count]
    coefficient_rank = length(coefficient_pivots)
    augmented_rank = length(all_pivots)
    consistent = coefficient_rank == augmented_rank
    pivot_set = Set(coefficient_pivots)
    free_columns = Int[j for j in 1:variable_count if !(j in pivot_set)]
    nullspace = _nullspace_from_rref(augmented, coefficient_pivots, variable_count)

    particular = if consistent
        value = zeros(ExactScalar, variable_count)
        for (pivot_row, pivot_column) in enumerate(coefficient_pivots)
            value[pivot_column] = augmented[pivot_row, variable_count + 1]
        end
        value
    else
        nothing
    end

    return (
        consistent=consistent,
        particular=particular,
        nullspace=nullspace,
        rank=coefficient_rank,
        augmented_rank=augmented_rank,
        pivot_columns=coefficient_pivots,
        free_columns=free_columns,
        rref=augmented,
    )
end

exact_linear_solve(
    A::AbstractMatrix,
    b::AbstractVector,
    policy::AbstractInexactPolicy,
) = exact_linear_solve(A, b; policy=policy)
