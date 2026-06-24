"""Abstract supertype for errors raised by RSDP."""
abstract type RSDPError <: Exception end

"""
    InexactDataError(value, policy, context, reason)

Raised when input data cannot be converted to an exact scalar under the
selected inexact-input policy.
"""
struct InexactDataError <: RSDPError
    value::Any
    policy::Any
    context::String
    reason::String
end

InexactDataError(value, policy, context::AbstractString, reason::AbstractString) =
    InexactDataError(value, policy, String(context), String(reason))

"""Compatibility alias for [`InexactDataError`](@ref)."""
const InexactConversionError = InexactDataError

function Base.showerror(io::IO, err::InexactDataError)
    print(
        io,
        "cannot exactify ",
        err.context,
        " = ",
        repr(err.value),
        " using ",
        nameof(typeof(err.policy)),
        ": ",
        err.reason,
    )
end

"""
    InvalidProblemError(message)

Raised when an exact conic problem has incompatible dimensions or malformed
data.
"""
struct InvalidProblemError <: RSDPError
    message::String
end

InvalidProblemError(message::AbstractString) = InvalidProblemError(String(message))

Base.showerror(io::IO, err::InvalidProblemError) =
    print(io, "invalid exact conic problem: ", err.message)

"""
    InconsistentAffineSystemError(rank, augmented_rank)

Raised when an API requiring an affine-space parameterization is called for an
inconsistent system.
"""
struct InconsistentAffineSystemError <: RSDPError
    rank::Int
    augmented_rank::Int
end

function Base.showerror(io::IO, err::InconsistentAffineSystemError)
    print(
        io,
        "affine system is inconsistent: rank(A) = ",
        err.rank,
        ", rank([A b]) = ",
        err.augmented_rank,
    )
end
