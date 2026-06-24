import SHA

export CERTIFICATE_VERSION,
    PrimalCertificate, ExactPrimalCertificate, make_primal_certificate

"""
Certificate format version emitted and accepted by this RSDP release.
"""
const CERTIFICATE_VERSION = v"0.1.0"

"""
    PrimalCertificate

An exact primal certificate tied to a particular conic problem.

`x` and the optional `objective` are always stored as `Rational{BigInt}`. A
certificate also records its format version, the exactification policy that
created the problem (when available), and untrusted descriptive metadata.

The problem hash is deliberately left in the representation supplied by the
problem implementation (for example, a hexadecimal string). Checkers never
trust the metadata.
"""
struct PrimalCertificate{H, P}
    certificate_version::VersionNumber
    problem_hash::H
    x::Vector{Rational{BigInt}}
    objective::Union{Nothing, Rational{BigInt}}
    exactification_policy::P
    metadata::Dict{Symbol, Any}
end

function PrimalCertificate(
    problem_hash,
    x::AbstractVector,
    objective = nothing;
    certificate_version::VersionNumber = CERTIFICATE_VERSION,
    exactification_policy = nothing,
    metadata::AbstractDict{Symbol} = Dict{Symbol, Any}(),
)
    exact_x = _certificate_exact_vector(x, "certificate primal vector")
    exact_objective =
        isnothing(objective) ? nothing :
        _certificate_exact_scalar(objective, "certificate objective")
    return PrimalCertificate{typeof(problem_hash), typeof(exactification_policy)}(
        certificate_version,
        problem_hash,
        exact_x,
        exact_objective,
        exactification_policy,
        Dict{Symbol, Any}(metadata),
    )
end

"""
Compatibility name for the package's exact primal certificate type.
"""
const ExactPrimalCertificate = PrimalCertificate

function Base.getproperty(certificate::PrimalCertificate, name::Symbol)
    if name === :x_exact || name === :exact_x
        return getfield(certificate, :x)
    elseif name === :objective_value
        return getfield(certificate, :objective)
    elseif name === :version
        return getfield(certificate, :certificate_version)
    end
    return getfield(certificate, name)
end

function Base.propertynames(::PrimalCertificate, private::Bool = false)
    names = (
        :certificate_version,
        :version,
        :problem_hash,
        :x,
        :x_exact,
        :objective,
        :objective_value,
        :exactification_policy,
        :metadata,
    )
    return private ? (names..., :exact_x) : names
end

struct _CertificateAutomaticObjective end
const _CERTIFICATE_AUTOMATIC_OBJECTIVE = _CertificateAutomaticObjective()
struct _CertificateAutomaticMetadata end
const _CERTIFICATE_AUTOMATIC_METADATA = _CertificateAutomaticMetadata()

"""
    make_primal_certificate(problem, x_exact; objective, include_objective=true)

Build a `PrimalCertificate` from an exact primal vector. Integer and rational
inputs are canonicalized to `Rational{BigInt}`.

By default the objective is computed exactly from the problem's objective
vector. Pass `include_objective=false` to omit that optional check, or pass an
explicit exact `objective` value to record it directly.

This constructor records evidence; use [`check_certificate`](@ref) to verify it.
"""
function make_primal_certificate(
    problem,
    x_exact::AbstractVector;
    objective = _CERTIFICATE_AUTOMATIC_OBJECTIVE,
    include_objective::Bool = true,
    metadata = _CERTIFICATE_AUTOMATIC_METADATA,
)
    x = _certificate_exact_vector(x_exact, "primal vector")
    problem_hash = _certificate_problem_hash(problem)

    exact_objective = if !include_objective
        nothing
    elseif objective isa _CertificateAutomaticObjective
        _certificate_objective_value(problem, x)
    elseif isnothing(objective)
        nothing
    else
        _certificate_exact_scalar(objective, "objective")
    end

    certificate_metadata = if metadata isa _CertificateAutomaticMetadata
        hasproperty(problem, :metadata) ?
        Dict{Symbol, Any}(getproperty(problem, :metadata)) : Dict{Symbol, Any}()
    elseif metadata isa AbstractDict{Symbol}
        Dict{Symbol, Any}(metadata)
    else
        throw(ArgumentError("certificate metadata must use Symbol keys"))
    end
    exactification_policy =
        hasproperty(problem, :exactification_policy) ?
        getproperty(problem, :exactification_policy) : nothing

    return PrimalCertificate(
        problem_hash,
        x,
        exact_objective;
        exactification_policy = exactification_policy,
        metadata = certificate_metadata,
    )
end

function _certificate_exact_scalar(value, label::AbstractString)
    if value isa Rational
        return BigInt(numerator(value)) // BigInt(denominator(value))
    elseif value isa Integer
        return BigInt(value) // BigInt(1)
    end
    throw(
        ArgumentError("$label must be an integer or rational value, got $(typeof(value))"),
    )
end

function _certificate_exact_vector(values::AbstractVector, label::AbstractString)
    result = Vector{Rational{BigInt}}(undef, length(values))
    for index in eachindex(values)
        result[index] = _certificate_exact_scalar(values[index], "$label entry $index")
    end
    return result
end

function _certificate_exact_matrix(values::AbstractMatrix, label::AbstractString)
    result = Matrix{Rational{BigInt}}(undef, size(values))
    for column in axes(values, 2), row in axes(values, 1)
        result[row, column] =
            _certificate_exact_scalar(values[row, column], "$label entry ($row, $column)")
    end
    return result
end

function _certificate_dot(coefficients::AbstractVector, x::AbstractVector)
    length(coefficients) == length(x) ||
        throw(DimensionMismatch("objective and primal vector lengths differ"))
    result = BigInt(0) // BigInt(1)
    for index in eachindex(coefficients, x)
        result +=
            _certificate_exact_scalar(coefficients[index], "objective coefficient") *
            _certificate_exact_scalar(x[index], "primal entry")
    end
    return result
end

function _certificate_problem_component(
    problem,
    names::Tuple{Vararg{Symbol}};
    required::Bool = true,
)
    for name in names
        if hasproperty(problem, name)
            return getproperty(problem, name)
        end
    end
    required || return nothing
    throw(
        ArgumentError(
            "problem $(typeof(problem)) does not expose any of $(collect(names))",
        ),
    )
end

_certificate_constraint_matrix(problem) =
    _certificate_problem_component(problem, (:A, :constraint_matrix))

_certificate_constraint_rhs(problem) =
    _certificate_problem_component(problem, (:b, :rhs, :constraint_rhs))

function _certificate_objective_coefficients(problem; required::Bool = true)
    return _certificate_problem_component(
        problem,
        (:c, :objective, :objective_coefficients);
        required = required,
    )
end

function _certificate_objective_constant(problem)
    if hasproperty(problem, :objective_constant)
        return _certificate_exact_scalar(
            getproperty(problem, :objective_constant),
            "objective constant",
        )
    elseif hasproperty(problem, :metadata)
        metadata = getproperty(problem, :metadata)
        if metadata isa AbstractDict && haskey(metadata, :objective_constant)
            return _certificate_exact_scalar(
                metadata[:objective_constant],
                "objective constant",
            )
        end
    end
    return zero(Rational{BigInt})
end

function _certificate_objective_value(problem, x)
    coefficients = _certificate_objective_coefficients(problem; required = false)
    isnothing(coefficients) && return nothing
    return _certificate_dot(coefficients, x) + _certificate_objective_constant(problem)
end

function _certificate_problem_hash(problem)
    # Integration seam: prefer the package's canonical hash function. Property
    # fallbacks support immutable problem records that cache their canonical
    # hash. The final fallback is deterministic and intended for lightweight
    # custom problem types used by downstream clients and tests.
    for function_name in (:problem_hash, :exact_problem_hash, :problem_fingerprint)
        if isdefined(@__MODULE__, function_name)
            hash_function = getfield(@__MODULE__, function_name)
            if hash_function isa Function && applicable(hash_function, problem)
                return hash_function(problem)
            end
        end
    end

    for property_name in (:problem_hash, :fingerprint, :hash)
        if hasproperty(problem, property_name)
            return getproperty(problem, property_name)
        end
    end

    return _certificate_stable_problem_fingerprint(problem)
end

function _certificate_stable_problem_fingerprint(problem)
    buffer = IOBuffer()
    _certificate_write_canonical(buffer, problem, IdDict{Any, Nothing}())
    bytes = take!(buffer)
    return bytes2hex(SHA.sha256(bytes))
end

function _certificate_write_canonical(io::IO, value, seen::IdDict{Any, Nothing})
    if value isa Rational
        print(io, "R", numerator(value), "/", denominator(value), ";")
    elseif value isa Integer
        print(io, "I", value, ";")
    elseif value isa AbstractFloat
        print(io, "F", bitstring(value), ";")
    elseif value isa Symbol
        print(io, "Y", String(value), ";")
    elseif value isa AbstractString
        print(io, "S", ncodeunits(value), ":", value, ";")
    elseif value isa Nothing
        print(io, "N;")
    elseif value isa AbstractDict
        entries = collect(pairs(value))
        sort!(entries; by = entry -> string(first(entry)))
        print(io, "D", string(typeof(value)), ":", length(entries), "{")
        for (key, entry_value) in entries
            _certificate_write_canonical(io, key, seen)
            _certificate_write_canonical(io, entry_value, seen)
        end
        print(io, "};")
    elseif value isa Tuple || value isa NamedTuple
        print(io, "C", string(typeof(value)), ":", length(value), "[")
        for entry in value
            _certificate_write_canonical(io, entry, seen)
        end
        print(io, "];")
    elseif value isa AbstractArray
        print(io, "C", string(typeof(value)), ":", size(value), "[")
        for entry in value
            _certificate_write_canonical(io, entry, seen)
        end
        print(io, "];")
    elseif isstructtype(typeof(value))
        if haskey(seen, value)
            throw(
                ArgumentError(
                    "cyclic problem data require a canonical problem_hash method",
                ),
            )
        end
        seen[value] = nothing
        print(io, "T", string(typeof(value)), "{")
        for name in fieldnames(typeof(value))
            field_name = name::Symbol
            print(io, string(field_name), "=")
            _certificate_write_canonical(io, getfield(value, field_name), seen)
        end
        print(io, "};")
        delete!(seen, value)
    else
        print(io, "O", string(typeof(value)), ":")
        show(io, value)
        print(io, ";")
    end
    return nothing
end
