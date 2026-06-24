"""
    AffineDiagnostics

Exact structural diagnostics for a linear system `A*x = b`.

`rank` is the rank of `A`; `augmented_rank` is the rank of `[A b]`.
Column indices are one-based and listed in increasing order.
"""
struct AffineDiagnostics
    row_count::Int
    variable_count::Int
    rank::Int
    augmented_rank::Int
    pivot_columns::Vector{Int}
    free_columns::Vector{Int}

    function AffineDiagnostics(
        row_count::Integer,
        variable_count::Integer,
        rank::Integer,
        augmented_rank::Integer,
        pivot_columns::AbstractVector{<:Integer},
        free_columns::AbstractVector{<:Integer},
    )
        m = Int(row_count)
        n = Int(variable_count)
        r = Int(rank)
        ar = Int(augmented_rank)
        m >= 0 || throw(ArgumentError("row_count must be nonnegative"))
        n >= 0 || throw(ArgumentError("variable_count must be nonnegative"))
        0 <= r <= min(m, n) || throw(ArgumentError("rank is out of range"))
        r <= ar <= min(m, n + 1) ||
            throw(ArgumentError("augmented_rank is out of range"))
        pivots = Int[p for p in pivot_columns]
        free = Int[p for p in free_columns]
        new(m, n, r, ar, pivots, free)
    end
end

"""Return the dimension of the nullspace recorded by `diagnostics`."""
nullity(diagnostics::AffineDiagnostics) =
    diagnostics.variable_count - diagnostics.rank

"""Return `true` when the recorded system is consistent."""
is_consistent(diagnostics::AffineDiagnostics) =
    diagnostics.rank == diagnostics.augmented_rank
