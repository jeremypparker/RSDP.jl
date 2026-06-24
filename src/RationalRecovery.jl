export NumericalPrimalHint,
    RecoveryOptions,
    RecoveryDiagnostics,
    RationalRecoveryError,
    recover_primal_certificate

"""
    NumericalPrimalHint(x; objective=nothing)

Floating-point primal information supplied to rational recovery. `objective` is
optional; when present, it is checked against the recovered exact objective.
"""
struct NumericalPrimalHint{T<:Real,O}
    x::Vector{T}
    objective::O
end

function Base.getproperty(hint::NumericalPrimalHint, name::Symbol)
    if name === :x_approx || name === :primal
        return getfield(hint, :x)
    elseif name === :objective_value
        return getfield(hint, :objective)
    end
    return getfield(hint, name)
end

function Base.propertynames(::NumericalPrimalHint, private::Bool = false)
    names = (:x, :objective, :x_approx, :primal, :objective_value)
    return names
end

function NumericalPrimalHint(x::AbstractVector{T}; objective = nothing) where {T<:Real}
    isempty(x) && return NumericalPrimalHint{T,typeof(objective)}(T[], objective)
    all(isfinite, x) ||
        throw(ArgumentError("numerical primal hint must contain only finite values"))
    if !isnothing(objective)
        objective isa Real ||
            throw(ArgumentError("objective hint must be real or nothing"))
        isfinite(objective) ||
            throw(ArgumentError("objective hint must be finite"))
    end
    return NumericalPrimalHint{T,typeof(objective)}(collect(x), objective)
end

NumericalPrimalHint(x::AbstractVector, objective) =
    NumericalPrimalHint(x; objective = objective)

"""
    RecoveryOptions(; max_denominator=1_000_000, atol=1e-8, rtol=1e-8,
                      objective_atol=atol, objective_rtol=rtol)

Controls rational recovery. `max_denominator` bounds the rationalized affine
coordinates `z`, not the entries of `x = xₚ + N*z`; the latter may inherit
denominators from the exact affine-space basis.
"""
struct RecoveryOptions
    max_denominator::BigInt
    atol::BigFloat
    rtol::BigFloat
    objective_atol::BigFloat
    objective_rtol::BigFloat
end

function RecoveryOptions(;
    max_denominator::Integer = 1_000_000,
    atol::Real = 1e-8,
    rtol::Real = 1e-8,
    objective_atol::Real = atol,
    objective_rtol::Real = rtol,
)
    max_denominator > 0 ||
        throw(ArgumentError("max_denominator must be positive"))
    for (name, value) in (
        (:atol, atol),
        (:rtol, rtol),
        (:objective_atol, objective_atol),
        (:objective_rtol, objective_rtol),
    )
        isfinite(value) && value >= 0 ||
            throw(ArgumentError("$name must be finite and nonnegative"))
    end
    return RecoveryOptions(
        BigInt(max_denominator),
        BigFloat(atol),
        BigFloat(rtol),
        BigFloat(objective_atol),
        BigFloat(objective_rtol),
    )
end

RecoveryOptions(max_denominator::Integer; kwargs...) =
    RecoveryOptions(; max_denominator = max_denominator, kwargs...)

"""
    RecoveryDiagnostics

Diagnostics for a recovery attempt. `stage` is `:success` on success and
identifies the failed boundary (`:input`, `:affine_space`, `:rationalization`,
`:proximity`, `:objective`, or `:certificate`) otherwise. `status` maps that
boundary to the package-wide [`ValidationStatus`](@ref).
"""
struct RecoveryDiagnostics
    stage::Symbol
    status::ValidationStatus
    message::String
    affine_dimension::Int
    max_coordinate_error::BigFloat
    max_primal_error::BigFloat
    certificate_report::Union{Nothing,CertificateCheckReport}
end

"""
    RationalRecoveryError

Exception raised when a numerical hint cannot be turned into a verified exact
certificate. Inspect its `diagnostics` field for the failed recovery boundary.
"""
struct RationalRecoveryError <: RSDPError
    diagnostics::RecoveryDiagnostics
end

function Base.showerror(io::IO, error::RationalRecoveryError)
    diagnostics = error.diagnostics
    print(
        io,
        "rational recovery failed at ",
        diagnostics.stage,
        ": ",
        diagnostics.message,
    )
    if isfinite(diagnostics.max_coordinate_error)
        print(io, " (max affine-coordinate error ", diagnostics.max_coordinate_error, ")")
    end
    if isfinite(diagnostics.max_primal_error)
        print(io, " (max primal error ", diagnostics.max_primal_error, ")")
    end
end

"""
    recover_primal_certificate(problem, hint; max_denominator=..., kwargs...)
    recover_primal_certificate(problem, hint, options; return_diagnostics=false)

Recover and verify an exact primal certificate from a numerical hint.

Recovery first obtains the exact affine space `x = xₚ + N*z`, estimates `z`
from the hint, rationalizes only those affine coordinates, and reconstructs
`x` exactly. This makes affine feasibility structural rather than a
floating-point afterthought. The candidate is then independently checked by
[`check_certificate`](@ref).

With `return_diagnostics=true`, return a named tuple containing both the
certificate and successful [`RecoveryDiagnostics`](@ref). Failures always
throw [`RationalRecoveryError`](@ref) with boundary-specific diagnostics.
"""
function recover_primal_certificate(
    problem,
    hint::NumericalPrimalHint;
    options::Union{Nothing,RecoveryOptions} = nothing,
    max_denominator::Integer = 1_000_000,
    atol::Real = 1e-8,
    rtol::Real = 1e-8,
    objective_atol::Real = atol,
    objective_rtol::Real = rtol,
    return_diagnostics::Bool = false,
)
    effective_options = isnothing(options) ? RecoveryOptions(;
        max_denominator = max_denominator,
        atol = atol,
        rtol = rtol,
        objective_atol = objective_atol,
        objective_rtol = objective_rtol,
    ) : options
    return recover_primal_certificate(
        problem,
        hint,
        effective_options;
        return_diagnostics = return_diagnostics,
    )
end

function recover_primal_certificate(
    problem,
    hint::NumericalPrimalHint,
    options::RecoveryOptions;
    return_diagnostics::Bool = false,
)
    x_particular, nullspace = try
        _recovery_exact_affine_space(problem)
    catch error
        _recovery_fail(
            :affine_space,
            "could not obtain exact affine space: $(sprint(showerror, error))",
        )
    end

    xp = try
        _certificate_exact_vector(x_particular, "affine particular solution")
    catch error
        _recovery_fail(:affine_space, sprint(showerror, error))
    end
    N = try
        _certificate_exact_matrix(nullspace, "affine nullspace basis")
    catch error
        _recovery_fail(:affine_space, sprint(showerror, error))
    end

    length(hint.x) == length(xp) ||
        _recovery_fail(
            :input,
            "hint has $(length(hint.x)) entries but affine space has $(length(xp))",
            size(N, 2),
        )
    size(N, 1) == length(xp) ||
        _recovery_fail(
            :affine_space,
            "nullspace basis has $(size(N, 1)) rows but xₚ has $(length(xp)) entries",
            size(N, 2),
        )

    z_approx = try
        _recovery_affine_coordinates(hint.x, xp, N)
    catch error
        _recovery_fail(
            :affine_space,
            "could not estimate affine coordinates: $(sprint(showerror, error))",
            size(N, 2),
        )
    end

    z_exact = Vector{Rational{BigInt}}(undef, length(z_approx))
    max_coordinate_error = BigFloat(0)
    for index in eachindex(z_approx)
        z_exact[index] = try
            _recovery_limit_denominator(z_approx[index], options.max_denominator)
        catch error
            _recovery_fail(
                :rationalization,
                "coordinate $index could not be rationalized: $(sprint(showerror, error))",
                length(z_approx),
                max_coordinate_error,
            )
        end
        error = abs(z_approx[index] - BigFloat(z_exact[index]))
        max_coordinate_error = max(max_coordinate_error, error)
    end

    x_exact = _recovery_reconstruct(xp, N, z_exact)
    max_primal_error, primal_scale = _recovery_primal_error(hint.x, x_exact)
    primal_tolerance = options.atol + options.rtol * primal_scale
    max_primal_error <= primal_tolerance ||
        _recovery_fail(
            :proximity,
            "nearest bounded-denominator affine point exceeds tolerance $primal_tolerance",
            length(z_exact),
            max_coordinate_error,
            max_primal_error,
        )

    certificate = make_primal_certificate(problem, x_exact)

    if !isnothing(hint.objective)
        if isnothing(certificate.objective)
            _recovery_fail(
                :objective,
                "hint contains an objective but the problem exposes no objective vector",
                length(z_exact),
                max_coordinate_error,
                max_primal_error,
            )
        end
        objective_error =
            abs(BigFloat(certificate.objective) - BigFloat(hint.objective))
        objective_scale = max(
            BigFloat(1),
            abs(BigFloat(certificate.objective)),
            abs(BigFloat(hint.objective)),
        )
        objective_tolerance =
            options.objective_atol + options.objective_rtol * objective_scale
        objective_error <= objective_tolerance ||
            _recovery_fail(
                :objective,
                "exact objective differs from hint by $objective_error " *
                "(tolerance $objective_tolerance)",
                length(z_exact),
                max_coordinate_error,
                max_primal_error,
            )
    end

    report = check_certificate(problem, certificate; diagnostics = true)
    report.valid ||
        _recovery_fail(
            :certificate,
            isempty(report.diagnostics) ?
            "independent certificate check failed" :
            join(report.diagnostics, "; "),
            length(z_exact),
            max_coordinate_error,
            max_primal_error,
            report,
        )

    diagnostics = RecoveryDiagnostics(
        :success,
        VALIDATED_PRIMAL_FEASIBLE,
        "recovered and independently verified an exact primal certificate",
        length(z_exact),
        max_coordinate_error,
        max_primal_error,
        report,
    )
    return return_diagnostics ?
           (certificate = certificate, diagnostics = diagnostics) :
           certificate
end

function _recovery_exact_affine_space(problem)
    space = if isdefined(@__MODULE__, :exact_affine_space)
        affine_function = getfield(@__MODULE__, :exact_affine_space)
        if applicable(affine_function, problem)
            affine_function(problem)
        elseif hasproperty(problem, :problem) &&
               applicable(affine_function, getproperty(problem, :problem))
            affine_function(getproperty(problem, :problem))
        else
            _recovery_solve_affine_fallback(problem)
        end
    else
        _recovery_solve_affine_fallback(problem)
    end

    if hasproperty(space, :solution)
        solution = getproperty(space, :solution)
        isnothing(solution) &&
            throw(
                ArgumentError(
                    "the problem's affine equations are not feasible",
                ),
            )
        space = solution
    end

    if space isa Tuple && length(space) == 2
        return space[1], space[2]
    end

    particular = nothing
    basis = nothing
    for name in (:xp, :x_p, :particular, :particular_solution, :x0)
        if hasproperty(space, name)
            particular = getproperty(space, name)
            break
        end
    end
    for name in (:N, :nullspace, :basis, :nullspace_basis)
        if hasproperty(space, name)
            basis = getproperty(space, name)
            break
        end
    end
    if isnothing(particular) || isnothing(basis)
        throw(
            ArgumentError(
                "exact_affine_space must return (xₚ, N) or expose particular " *
                "and nullspace fields",
            ),
        )
    end
    return particular, basis
end

function _recovery_solve_affine_fallback(problem)
    isdefined(@__MODULE__, :solve_affine) ||
        throw(ArgumentError("the package must provide exact_affine_space"))
    affine_function = getfield(@__MODULE__, :solve_affine)
    if applicable(affine_function, problem)
        return affine_function(problem)
    elseif hasproperty(problem, :problem) &&
           applicable(affine_function, getproperty(problem, :problem))
        return affine_function(getproperty(problem, :problem))
    end
    throw(MethodError(affine_function, (problem,)))
end

function _recovery_affine_coordinates(hint, xp, N)
    rows, columns = size(N)
    columns == 0 && return BigFloat[]

    gram = zeros(BigFloat, columns, columns)
    rhs = zeros(BigFloat, columns)
    for row in 1:rows
        displacement = BigFloat(hint[row]) - BigFloat(xp[row])
        for left in 1:columns
            n_left = BigFloat(N[row, left])
            rhs[left] += n_left * displacement
            for right in left:columns
                gram[left, right] += n_left * BigFloat(N[row, right])
            end
        end
    end
    for left in 1:columns, right in 1:(left - 1)
        gram[left, right] = gram[right, left]
    end
    return _recovery_solve_dense(gram, rhs)
end

function _recovery_solve_dense(matrix::Matrix{BigFloat}, rhs::Vector{BigFloat})
    n = length(rhs)
    size(matrix) == (n, n) ||
        throw(DimensionMismatch("linear system must be square"))
    augmented = hcat(copy(matrix), copy(rhs))

    for pivot_column in 1:n
        pivot_row = pivot_column
        pivot_size = abs(augmented[pivot_row, pivot_column])
        for row in (pivot_column + 1):n
            candidate_size = abs(augmented[row, pivot_column])
            if candidate_size > pivot_size
                pivot_row = row
                pivot_size = candidate_size
            end
        end
        iszero(pivot_size) &&
            throw(ArgumentError("nullspace basis columns are linearly dependent"))
        if pivot_row != pivot_column
            for column in pivot_column:(n + 1)
                augmented[pivot_column, column], augmented[pivot_row, column] =
                    augmented[pivot_row, column], augmented[pivot_column, column]
            end
        end

        pivot = augmented[pivot_column, pivot_column]
        for row in (pivot_column + 1):n
            factor = augmented[row, pivot_column] / pivot
            augmented[row, pivot_column] = 0
            for column in (pivot_column + 1):(n + 1)
                augmented[row, column] -= factor * augmented[pivot_column, column]
            end
        end
    end

    solution = zeros(BigFloat, n)
    for row in n:-1:1
        remainder = augmented[row, n + 1]
        for column in (row + 1):n
            remainder -= augmented[row, column] * solution[column]
        end
        solution[row] = remainder / augmented[row, row]
    end
    return solution
end

function _recovery_limit_denominator(value::Real, max_denominator::BigInt)
    isfinite(value) || throw(ArgumentError("value must be finite"))
    exact = rationalize(BigInt, BigFloat(value); tol = BigFloat(0))
    denominator(exact) <= max_denominator && return exact

    sign = exact < 0 ? -BigInt(1) : BigInt(1)
    numerator_remaining = abs(numerator(exact))
    denominator_remaining = denominator(exact)
    p0, q0 = BigInt(0), BigInt(1)
    p1, q1 = BigInt(1), BigInt(0)

    while true
        quotient = div(numerator_remaining, denominator_remaining)
        q2 = q0 + quotient * q1
        q2 > max_denominator && break
        p0, q0, p1, q1 =
            p1, q1, p0 + quotient * p1, q2
        numerator_remaining, denominator_remaining =
            denominator_remaining,
            numerator_remaining - quotient * denominator_remaining
    end

    multiplier = div(max_denominator - q0, q1)
    bound1 = (p0 + multiplier * p1) // (q0 + multiplier * q1)
    bound2 = p1 // q1
    positive_exact = abs(exact)
    nearest =
        abs(bound2 - positive_exact) <= abs(bound1 - positive_exact) ?
        bound2 :
        bound1
    return sign * nearest
end

function _recovery_reconstruct(xp, N, z)
    result = copy(xp)
    for row in axes(N, 1), column in axes(N, 2)
        result[row] += N[row, column] * z[column]
    end
    return result
end

function _recovery_primal_error(hint, exact)
    maximum_error = BigFloat(0)
    scale = BigFloat(1)
    for index in eachindex(hint, exact)
        numerical = BigFloat(hint[index])
        exact_value = BigFloat(exact[index])
        maximum_error = max(maximum_error, abs(numerical - exact_value))
        scale = max(scale, abs(numerical), abs(exact_value))
    end
    return maximum_error, scale
end

function _recovery_fail(
    stage::Symbol,
    message::AbstractString,
    affine_dimension::Integer = 0,
    max_coordinate_error::Real = BigFloat(Inf),
    max_primal_error::Real = BigFloat(Inf),
    report::Union{Nothing,CertificateCheckReport} = nothing,
)
    throw(
        RationalRecoveryError(
            RecoveryDiagnostics(
                stage,
                _recovery_failure_status(stage, report),
                String(message),
                Int(affine_dimension),
                BigFloat(max_coordinate_error),
                BigFloat(max_primal_error),
                report,
            ),
        ),
    )
end

function _recovery_failure_status(
    stage::Symbol,
    report::Union{Nothing,CertificateCheckReport},
)
    if stage === :affine_space
        return RECOVERY_FAILED_AFFINE
    elseif stage === :rationalization || stage === :proximity
        return RECOVERY_FAILED_DENOMINATOR_LIMIT
    elseif stage === :certificate
        if !isnothing(report) && !report.cones_valid
            return RECOVERY_FAILED_CONE
        end
        return CERTIFICATE_CHECK_FAILED
    end
    return RECOVERY_FAILED
end
