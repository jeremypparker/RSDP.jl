"""Exact scalar type used throughout RSDP."""
const ExactScalar = Rational{BigInt}

"""Abstract supertype for policies controlling conversion of inexact data."""
abstract type AbstractInexactPolicy end

"""Compatibility alias for [`AbstractInexactPolicy`](@ref)."""
const InexactPolicy = AbstractInexactPolicy

"""
    ErrorOnInexact()

Reject every floating-point or otherwise inexact scalar. This is the default
policy for all public exactification and problem-construction entry points.
"""
struct ErrorOnInexact <: AbstractInexactPolicy end

"""
    RationalizeInexact(; tolerance=nothing)
    RationalizeInexact(tolerance)

Convert finite floating-point values with Julia's continued-fraction
`rationalize`. When `tolerance` is omitted, Julia's deterministic
type-appropriate default tolerance is used.
"""
struct RationalizeInexact{T} <: AbstractInexactPolicy
    tolerance::T
end

function RationalizeInexact(; tolerance=nothing, tol=nothing)
    if !isnothing(tolerance) && !isnothing(tol)
        throw(ArgumentError("specify only one of tolerance and tol"))
    end
    selected_tolerance = isnothing(tolerance) ? tol : tolerance
    isnothing(selected_tolerance) && return RationalizeInexact{Nothing}(nothing)
    selected_tolerance isa Real ||
        throw(ArgumentError("rationalization tolerance must be real"))
    return RationalizeInexact(selected_tolerance)
end

function RationalizeInexact(tolerance::Real)
    isfinite(tolerance) ||
        throw(ArgumentError("rationalization tolerance must be finite"))
    tolerance >= zero(tolerance) ||
        throw(ArgumentError("rationalization tolerance must be nonnegative"))
    return RationalizeInexact{typeof(tolerance)}(tolerance)
end

"""
    DecimalStringInexact()

Convert a finite floating-point value by parsing its shortest decimal string.
For example, `0.1` becomes exactly `1//10`, rather than the exact rational
representation of its binary floating-point payload.
"""
struct DecimalStringInexact <: AbstractInexactPolicy end

"""The strict policy used when no inexact-input policy is supplied."""
const DEFAULT_INEXACT_POLICY = ErrorOnInexact()

"""
    exactify(value[, policy]; context="value")

Convert a scalar or array to [`ExactScalar`](@ref). Integer, rational, and
decimal-string inputs are exact and are always accepted. Inexact numeric input
is handled according to `policy`, which defaults to [`ErrorOnInexact`](@ref).
"""
exactify(value; policy::AbstractInexactPolicy=DEFAULT_INEXACT_POLICY, context="value") =
    exactify(value, policy; context=context)

exactify(value::ExactScalar, ::AbstractInexactPolicy; context="value") = value

exactify(value::Integer, ::AbstractInexactPolicy; context="value") =
    BigInt(value) // BigInt(1)

exactify(value::Rational, ::AbstractInexactPolicy; context="value") =
    BigInt(numerator(value)) // BigInt(denominator(value))

function exactify(
    value::AbstractFloat,
    policy::ErrorOnInexact;
    context="value",
)
    throw(
        InexactDataError(
            value,
            policy,
            context,
            "floating-point input is forbidden by the strict default policy",
        ),
    )
end

function exactify(
    value::AbstractFloat,
    policy::RationalizeInexact;
    context="value",
)
    isfinite(value) ||
        throw(InexactDataError(value, policy, context, "value must be finite"))
    value == zero(value) && return zero(ExactScalar)
    if isnothing(policy.tolerance)
        return rationalize(BigInt, value)
    end
    return rationalize(BigInt, value; tol=policy.tolerance)
end

function exactify(
    value::AbstractFloat,
    policy::DecimalStringInexact;
    context="value",
)
    isfinite(value) ||
        throw(InexactDataError(value, policy, context, "value must be finite"))
    return _parse_exact_decimal(string(value), policy, context)
end

function exactify(
    value::AbstractString,
    policy::AbstractInexactPolicy;
    context="value",
)
    return _parse_exact_decimal(value, policy, context)
end

function exactify(value::Real, policy::AbstractInexactPolicy; context="value")
    throw(
        InexactDataError(
            value,
            policy,
            context,
            "unsupported real scalar type $(typeof(value))",
        ),
    )
end

function exactify(value, policy::AbstractInexactPolicy; context="value")
    throw(
        InexactDataError(
            value,
            policy,
            context,
            "expected an integer, rational, decimal string, or finite float",
        ),
    )
end

function exactify(
    values::AbstractArray,
    policy::AbstractInexactPolicy;
    context="array",
)
    result = Array{ExactScalar}(undef, size(values))
    for index in CartesianIndices(values)
        location = string(context, "[", join(Tuple(index), ","), "]")
        result[index] = exactify(values[index], policy; context=location)
    end
    return result
end

function _parse_exact_decimal(
    source::AbstractString,
    policy::AbstractInexactPolicy,
    context,
)
    text = strip(source)
    isempty(text) &&
        throw(InexactDataError(source, policy, context, "empty numeric string"))

    slash_parts = split(text, '/'; limit=2)
    if length(slash_parts) == 2
        try
            numerator_value = parse(BigInt, strip(slash_parts[1]))
            denominator_value = parse(BigInt, strip(slash_parts[2]))
            denominator_value == 0 &&
                throw(
                    InexactDataError(
                        source,
                        policy,
                        context,
                        "rational denominator must be nonzero",
                    ),
                )
            return numerator_value // denominator_value
        catch err
            err isa InexactDataError && rethrow()
            throw(
                InexactDataError(
                    source,
                    policy,
                    context,
                    "invalid rational string",
                ),
            )
        end
    end

    matched = match(
        r"^([+-]?)(?:(\d+)(?:\.(\d*))?|\.(\d+))(?:[eE]([+-]?\d+))?$",
        text,
    )
    isnothing(matched) &&
        throw(InexactDataError(source, policy, context, "invalid decimal string"))

    sign_text, integer_digits, fractional_digits_1, fractional_digits_2, exponent_text =
        matched.captures
    integer_part = isnothing(integer_digits) ? "0" : integer_digits
    fractional_part = if !isnothing(fractional_digits_2)
        fractional_digits_2
    elseif isnothing(fractional_digits_1)
        ""
    else
        fractional_digits_1
    end
    exponent = isnothing(exponent_text) ? 0 : try
        parse(Int, exponent_text)
    catch
        throw(InexactDataError(source, policy, context, "decimal exponent is out of range"))
    end

    coefficient = parse(BigInt, string(integer_part, fractional_part))
    sign_text == "-" && (coefficient = -coefficient)
    decimal_places = length(fractional_part) - exponent
    if decimal_places >= 0
        return coefficient // (BigInt(10)^decimal_places)
    end
    return (coefficient * BigInt(10)^(-decimal_places)) // BigInt(1)
end
