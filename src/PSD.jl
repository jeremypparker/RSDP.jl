if !isdefined(@__MODULE__, :ExactRational)
    const ExactRational = Rational{BigInt}
end

export PSDCheckDiagnostic,
    PSDCheckResult,
    check_psd,
    check_psd_exact,
    exact_psd_check,
    is_positive_semidefinite_exact,
    is_psd_exact,
    is_psd

"""
    PSDCheckDiagnostic

A machine-readable diagnostic from [`check_psd`](@ref). `code` identifies the
outcome, while `step` and `index` locate an elimination pivot when applicable.
"""
struct PSDCheckDiagnostic
    code::Symbol
    message::String
    step::Int
    index::Union{Nothing, Int}
end

"""
    PSDCheckResult

Structured result of an exact positive-semidefinite check.

Fields:

  - `is_psd`: whether the input is exactly symmetric positive semidefinite;
  - `rank`: exact rank accumulated before the outcome was established;
  - `diagnostic`: machine-readable status or failure information;
  - `permutation`: pivot order, followed by uneliminated input indices;
  - `pivots`: exact positive pivots, followed by zero pivots on success;
  - `witness`: an exact vector `x` with `x' * A * x < 0` for algebraic
    indefiniteness, when available.
"""
struct PSDCheckResult
    is_psd::Bool
    rank::Int
    diagnostic::PSDCheckDiagnostic
    permutation::Vector{Int}
    pivots::Vector{ExactRational}
    witness::Union{Nothing, Vector{ExactRational}}
end

function Base.show(io::IO, result::PSDCheckResult)
    print(
        io,
        "PSDCheckResult(",
        result.is_psd ? "positive semidefinite" : "not positive semidefinite",
        ", rank=",
        result.rank,
        ", code=:",
        result.diagnostic.code,
        ")",
    )
end

_psd_exact_rational(value::ExactRational) = value
_psd_exact_rational(value::Integer) = BigInt(value) // BigInt(1)
_psd_exact_rational(value::Rational{T}) where {T <: Integer} =
    BigInt(numerator(value)) // BigInt(denominator(value))
_psd_exact_rational(value) = nothing

function _exact_matrix(matrix::AbstractMatrix)
    rows, columns = size(matrix)
    exact = Matrix{ExactRational}(undef, rows, columns)
    for column in 1:columns
        for row in 1:rows
            value = _psd_exact_rational(matrix[row, column])
            value === nothing && return nothing, (row, column), typeof(matrix[row, column])
            exact[row, column] = value
        end
    end
    return exact, nothing, nothing
end

function _psd_result(
    is_psd::Bool,
    rank::Int,
    code::Symbol,
    message::String,
    step::Int,
    index::Union{Nothing, Int},
    permutation::Vector{Int},
    pivots::Vector{ExactRational},
    witness::Union{Nothing, Vector{ExactRational}} = nothing,
)
    diagnostic = PSDCheckDiagnostic(code, message, step, index)
    return PSDCheckResult(is_psd, rank, diagnostic, permutation, pivots, witness)
end

"""
    check_psd(matrix) -> PSDCheckResult

Check positive semidefiniteness using exact rational arithmetic only.

The algorithm performs symmetric positive-pivot elimination. At each step,
`A ⪰ 0` is equivalent to a positive pivot together with its exact Schur
complement being PSD. If no positive diagonal remains, PSD requires the entire
remainder to be zero. This criterion handles singular matrices and zero leading
pivots without division by zero or floating-point eigenvalue tolerances.

Integer and rational matrices are accepted and converted to
`Rational{BigInt}`. Floating-point or otherwise inexact entries, nonsquare
matrices, and nonsymmetric matrices produce structured negative results.
"""
function check_psd(matrix::AbstractMatrix)
    rows, columns = size(matrix)
    if rows != columns
        return _psd_result(
            false,
            0,
            :nonsquare,
            "PSD checking requires a square matrix, got size $(size(matrix))",
            0,
            nothing,
            Int[],
            ExactRational[],
        )
    end

    exact, bad_index, bad_type = _exact_matrix(matrix)
    if exact === nothing
        row, column = bad_index
        return _psd_result(
            false,
            0,
            :inexact_input,
            "entry ($row, $column) has non-exact type $bad_type; use integers or rationals",
            0,
            row,
            collect(1:rows),
            ExactRational[],
        )
    end

    for column in 1:columns
        for row in 1:(column-1)
            if exact[row, column] != exact[column, row]
                return _psd_result(
                    false,
                    0,
                    :nonsymmetric,
                    "entries ($row, $column) and ($column, $row) are unequal",
                    0,
                    row,
                    collect(1:rows),
                    ExactRational[],
                )
            end
        end
    end

    n = rows
    remainder = copy(exact)
    basis = zeros(ExactRational, n, n)
    for index in 1:n
        basis[index, index] = one(ExactRational)
    end
    labels = collect(1:n)
    pivot_order = Int[]
    pivots = ExactRational[]
    rank = 0

    while !isempty(labels)
        remainder_size = length(labels)

        negative_index = findfirst(i -> remainder[i, i] < 0, 1:remainder_size)
        if negative_index !== nothing
            original_index = labels[negative_index]
            witness = copy(basis[:, negative_index])
            permutation = vcat(pivot_order, labels)
            return _psd_result(
                false,
                rank,
                :negative_pivot,
                "exact congruence reduction found a negative diagonal entry",
                rank + 1,
                original_index,
                permutation,
                copy(pivots),
                witness,
            )
        end

        positive_index = findfirst(i -> remainder[i, i] > 0, 1:remainder_size)
        if positive_index === nothing
            off_diagonal = nothing
            for column in 2:remainder_size
                for row in 1:(column-1)
                    if !iszero(remainder[row, column])
                        off_diagonal = (row, column)
                        break
                    end
                end
                off_diagonal === nothing || break
            end

            if off_diagonal !== nothing
                row, column = off_diagonal
                value = remainder[row, column]
                sign = value > 0 ? one(ExactRational) : -one(ExactRational)
                witness = basis[:, row] - sign * basis[:, column]
                permutation = vcat(pivot_order, labels)
                return _psd_result(
                    false,
                    rank,
                    :zero_diagonal_nonzero_row,
                    "a zero-diagonal remainder has a nonzero off-diagonal entry",
                    rank + 1,
                    labels[row],
                    permutation,
                    copy(pivots),
                    witness,
                )
            end

            append!(pivot_order, labels)
            append!(pivots, fill(zero(ExactRational), remainder_size))
            return _psd_result(
                true,
                rank,
                rank == n ? :positive_definite : :positive_semidefinite,
                rank == n ? "all exact pivots are positive" :
                "the exact Schur remainder is zero",
                rank,
                nothing,
                pivot_order,
                pivots,
            )
        end

        pivot = remainder[positive_index, positive_index]
        push!(pivot_order, labels[positive_index])
        push!(pivots, pivot)
        rank += 1

        remaining = [i for i in 1:remainder_size if i != positive_index]
        if isempty(remaining)
            empty!(labels)
            remainder = Matrix{ExactRational}(undef, 0, 0)
            basis = Matrix{ExactRational}(undef, n, 0)
            continue
        end

        coupling = copy(remainder[remaining, positive_index])
        remainder =
            remainder[remaining, remaining] - (coupling * transpose(coupling)) / pivot
        pivot_basis = copy(basis[:, positive_index])
        basis = basis[:, remaining] - pivot_basis * transpose(coupling / pivot)
        labels = labels[remaining]
    end

    return _psd_result(
        true,
        rank,
        :positive_definite,
        "all exact pivots are positive",
        rank,
        nothing,
        pivot_order,
        pivots,
    )
end

"""
Alias for [`check_psd`](@ref).
"""
exact_psd_check(matrix::AbstractMatrix) = check_psd(matrix)

"""
    is_positive_semidefinite_exact(matrix)

Return only the Boolean outcome of [`check_psd`](@ref).
"""
is_positive_semidefinite_exact(matrix::AbstractMatrix) = check_psd(matrix).is_psd

"""
Short Boolean alias for [`is_positive_semidefinite_exact`](@ref).
"""
is_psd_exact(matrix::AbstractMatrix) = is_positive_semidefinite_exact(matrix)

"""
Compatibility alias for [`check_psd`](@ref).
"""
check_psd_exact(matrix::AbstractMatrix) = check_psd(matrix)

"""
Conventional Boolean alias for [`is_positive_semidefinite_exact`](@ref).
"""
is_psd(matrix::AbstractMatrix) = is_positive_semidefinite_exact(matrix)
