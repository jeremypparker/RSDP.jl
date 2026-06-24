export CertificateCheckReport, ValidationReport, check_certificate

"""
    CertificateCheckReport

Detailed result of independently checking a primal certificate. The individual
flags cover the problem fingerprint, affine equations, cone membership, and the
optional objective value.
"""
struct CertificateCheckReport
    valid::Bool
    status::ValidationStatus
    version_valid::Bool
    problem_hash_valid::Bool
    affine_valid::Bool
    cones_valid::Bool
    objective_valid::Bool
    diagnostics::Vector{String}
    computed_objective::Union{Nothing, Rational{BigInt}}
end

"""
Compatibility name for the authoritative certificate validation report.
"""
const ValidationReport = CertificateCheckReport

function Base.getproperty(report::CertificateCheckReport, name::Symbol)
    if name === :ok
        return getfield(report, :valid)
    end
    return getfield(report, name)
end

function Base.propertynames(::CertificateCheckReport, private::Bool = false)
    return (
        :ok,
        :valid,
        :status,
        :version_valid,
        :problem_hash_valid,
        :affine_valid,
        :cones_valid,
        :objective_valid,
        :diagnostics,
        :computed_objective,
    )
end

"""
    check_certificate(problem, certificate; diagnostics=false, check_objective=true)

Independently verify a `PrimalCertificate` using exact arithmetic.

The checker recomputes the problem hash, verifies `A * x == b`, delegates exact
block membership to the package's `check_cones` API, and, when present, checks
the recorded objective. It never trusts feasibility metadata from the
certificate itself.

The authoritative [`ValidationReport`](@ref) is returned by default. Set
`diagnostics=false` when only its Boolean outcome is needed.
"""
function check_certificate(
    problem,
    certificate::PrimalCertificate;
    diagnostics::Bool = true,
    check_objective::Bool = true,
)
    messages = String[]
    version_valid = certificate.certificate_version == CERTIFICATE_VERSION
    version_valid || push!(
        messages,
        "unsupported certificate version $(certificate.certificate_version); " *
        "expected $CERTIFICATE_VERSION",
    )

    current_hash = _certificate_problem_hash(problem)
    hash_valid = isequal(current_hash, certificate.problem_hash)
    hash_valid || push!(
        messages,
        "problem hash mismatch: certificate is for a different or modified problem",
    )

    affine_valid = _certificate_check_affine(problem, certificate.x, messages)
    cones_valid = _certificate_check_cones(problem, certificate.x, messages)

    computed_objective = nothing
    objective_valid = true
    if check_objective && !isnothing(certificate.objective)
        coefficients = _certificate_objective_coefficients(problem; required = false)
        if isnothing(coefficients)
            objective_valid = false
            push!(
                messages,
                "certificate records an objective but the problem exposes " *
                "no objective vector",
            )
        else
            try
                computed_objective = _certificate_objective_value(problem, certificate.x)
                objective_valid = computed_objective == certificate.objective
                objective_valid || push!(
                    messages,
                    "objective mismatch: expected $(certificate.objective), " *
                    "computed $computed_objective",
                )
            catch error
                objective_valid = false
                push!(messages, "objective check failed: $(sprint(showerror, error))")
            end
        end
    end

    valid = version_valid && hash_valid && affine_valid && cones_valid && objective_valid
    status = valid ? VALIDATED_PRIMAL_FEASIBLE : CERTIFICATE_CHECK_FAILED
    report = CertificateCheckReport(
        valid,
        status,
        version_valid,
        hash_valid,
        affine_valid,
        cones_valid,
        objective_valid,
        messages,
        computed_objective,
    )
    return diagnostics ? report : report.valid
end

function _certificate_check_affine(
    problem,
    x::Vector{Rational{BigInt}},
    messages::Vector{String},
)
    try
        A = _certificate_constraint_matrix(problem)
        b = _certificate_constraint_rhs(problem)
        size(A, 2) == length(x) || begin
            push!(
                messages,
                "affine dimension mismatch: A has $(size(A, 2)) columns but " *
                "x has $(length(x)) entries",
            )
            return false
        end
        size(A, 1) == length(b) || begin
            push!(
                messages,
                "affine dimension mismatch: A has $(size(A, 1)) rows but " *
                "b has $(length(b)) entries",
            )
            return false
        end

        for row in axes(A, 1)
            residual = -_certificate_exact_scalar(b[row], "constraint rhs")
            for column in axes(A, 2)
                residual +=
                    _certificate_exact_scalar(A[row, column], "constraint coefficient") *
                    x[column]
            end
            if !iszero(residual)
                push!(messages, "affine equation $row failed with exact residual $residual")
                return false
            end
        end
        return true
    catch error
        push!(messages, "affine check failed: $(sprint(showerror, error))")
        return false
    end
end

function _certificate_check_cones(
    problem,
    x::Vector{Rational{BigInt}},
    messages::Vector{String},
)
    try
        result = _certificate_call_check_cones(problem, x)
        valid, detail = _certificate_cone_result(result)
        if !valid
            suffix = isempty(detail) ? "" : ": $detail"
            push!(messages, "cone membership check failed$suffix")
        end
        return valid
    catch error
        push!(messages, "cone membership check failed: $(sprint(showerror, error))")
        return false
    end
end

function _certificate_call_check_cones(problem, x)
    cones = _certificate_problem_component(
        problem,
        (:cones, :cone, :cone_blocks, :blocks, :K);
        required = false,
    )
    if hasproperty(problem, :cones) && isnothing(getproperty(problem, :cones))
        return true
    end

    if isdefined(@__MODULE__, :check_cones)
        cone_checker = getfield(@__MODULE__, :check_cones)

        # Primary integration contract.
        applicable(cone_checker, problem, x) && return cone_checker(problem, x)

        # Narrow compatibility fallbacks for a checker organized around cone blocks.
        if !isnothing(cones)
            applicable(cone_checker, cones, x) && return cone_checker(cones, x)
            applicable(cone_checker, x, cones) && return cone_checker(x, cones)
        end
    end

    # The current low-level cone API is block-oriented. Supporting it here keeps
    # the checker usable for problem wrappers while `check_cones(problem, x)` is
    # the preferred package-level integration point.
    if !isnothing(cones) && isdefined(@__MODULE__, :check_cone_membership)
        block_checker = getfield(@__MODULE__, :check_cone_membership)
        applicable(block_checker, cones, x) && return block_checker(cones, x)
    end

    throw(
        ArgumentError(
            "the package must provide check_cones(problem, x), or the problem " *
            "must expose cone blocks accepted by check_cone_membership",
        ),
    )
end

function _certificate_cone_result(result::ConeMembershipResult)
    detail = join((diagnostic.message for diagnostic in result.diagnostics), "; ")
    return result.is_member, detail
end

function _certificate_cone_result(result)
    if result isa Bool
        return result, ""
    elseif result isa AbstractVector
        pairs = map(_certificate_cone_result, result)
        return all(first, pairs), join(filter(x -> !isempty(x), last.(pairs)), "; ")
    elseif result isa Tuple
        isempty(result) && return true, ""
        if first(result) isa Bool
            detail = length(result) >= 2 ? string(result[2]) : ""
            return first(result), detail
        end
        pairs = map(_certificate_cone_result, result)
        return all(first, pairs), join(filter(x -> !isempty(x), last.(pairs)), "; ")
    end

    for name in (:valid, :ok, :feasible, :success, :is_member)
        if hasproperty(result, name)
            valid = Bool(getproperty(result, name))
            detail = if hasproperty(result, :diagnostics)
                string(getproperty(result, :diagnostics))
            elseif hasproperty(result, :message)
                string(getproperty(result, :message))
            else
                ""
            end
            return valid, detail
        end
    end

    throw(
        ArgumentError(
            "unsupported check_cones result $(typeof(result)); expected Bool " *
            "or an object with a validity field",
        ),
    )
end
