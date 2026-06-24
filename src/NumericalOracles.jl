"""
    AbstractNumericalOracle

Supertype for numerical solvers that produce non-authoritative primal hints.
Exact validation remains the responsibility of rational recovery and
[`check_certificate`](@ref).
"""
abstract type AbstractNumericalOracle end

"""
    HypatiaOracle(; optimizer_factory=nothing, attributes=Dict(), silent=true)

Configuration for the optional Hypatia numerical oracle. Hypatia is a weak
dependency: `using Hypatia` loads the implementation in `RSDPHypatiaExt`.
The optimizer only generates a floating-point hint and is never called by the
exact certificate checker.
"""
struct HypatiaOracle <: AbstractNumericalOracle
    optimizer_factory::Any
    attributes::Dict{Any, Any}
    silent::Bool
end

function HypatiaOracle(;
    optimizer_factory = nothing,
    attributes = Dict{Any, Any}(),
    silent::Bool = true,
)
    return HypatiaOracle(optimizer_factory, Dict{Any, Any}(attributes), silent)
end

"""
    NumericalOracleResult

Result of numerical hint generation. `status` is a package-wide validation
status, `primal` and `objective_value` are floating-point hints when available,
and `raw_status`/`raw_result` preserve solver information for diagnostics.
No field in this type is an exact mathematical claim.
"""
struct NumericalOracleResult{T}
    status::ValidationStatus
    primal::Union{Nothing, Vector{T}}
    objective_value::Union{Nothing, T}
    raw_status::Any
    raw_result::Any
    diagnostics::Vector{String}
end

function NumericalOracleResult(
    status::ValidationStatus;
    primal = nothing,
    objective_value = nothing,
    raw_status = nothing,
    raw_result = nothing,
    diagnostics = String[],
)
    T = if !isnothing(primal)
        eltype(primal)
    elseif !isnothing(objective_value)
        typeof(objective_value)
    else
        Float64
    end
    converted_primal = isnothing(primal) ? nothing : Vector{T}(primal)
    converted_objective = isnothing(objective_value) ? nothing : convert(T, objective_value)
    return NumericalOracleResult{T}(
        status,
        converted_primal,
        converted_objective,
        raw_status,
        raw_result,
        String[String(message) for message in diagnostics],
    )
end

"""
    OracleValidationResult

Combined numerical-oracle, rational-recovery, and exact-checking result.
Only `report.ok == true` represents validated primal feasibility.
"""
struct OracleValidationResult{T}
    oracle_result::NumericalOracleResult{T}
    certificate::Union{Nothing, ExactPrimalCertificate}
    report::ValidationReport
    diagnostics::Vector{String}
end

"""
    solve_oracle(problem, oracle; kwargs...)

Generate a numerical primal hint for `problem`. Concrete oracle packages add
methods for this function. The core fallback returns a structured failure.
"""
function solve_oracle(
    problem::ExactConicProblem,
    oracle::AbstractNumericalOracle;
    kwargs...,
)
    message = "no numerical implementation is loaded for oracle $(typeof(oracle))"
    return NumericalOracleResult(
        NUMERICAL_ORACLE_FAILED;
        raw_status = :implementation_not_loaded,
        diagnostics = [message],
    )
end

function _failed_validation_report(status::ValidationStatus, diagnostics::Vector{String})
    return CertificateCheckReport(
        false,
        status,
        false,
        false,
        false,
        false,
        false,
        copy(diagnostics),
        nothing,
    )
end

"""
    validate_with_oracle(problem, oracle; kwargs...)

Ask `oracle` for a floating-point primal hint, recover an exact rational
certificate using [`recover_primal_certificate`](@ref), and independently
check that certificate. Failures are returned as structured diagnostics rather
than promoted to exact claims.

The recovery keywords match [`RecoveryOptions`](@ref). This function validates
primal feasibility only; it does not claim optimality.
"""
function validate_with_oracle(
    problem::ExactConicProblem,
    oracle::AbstractNumericalOracle;
    options::Union{Nothing,RecoveryOptions} = nothing,
    max_denominator::Integer = 1_000_000,
    atol::Real = 1e-8,
    rtol::Real = 1e-8,
    objective_atol::Real = atol,
    objective_rtol::Real = rtol,
)
    oracle_result = try
        solve_oracle(problem, oracle)
    catch error
        message = "numerical oracle threw an exception: $(sprint(showerror, error))"
        NumericalOracleResult(
            NUMERICAL_ORACLE_FAILED;
            raw_status = :exception,
            raw_result = error,
            diagnostics = [message],
        )
    end
    diagnostics = copy(oracle_result.diagnostics)

    if isnothing(oracle_result.primal)
        isempty(diagnostics) &&
            push!(diagnostics, "numerical oracle returned no primal vector")
        report = _failed_validation_report(
            oracle_result.status == UNSUPPORTED_CONE ? UNSUPPORTED_CONE :
            NUMERICAL_ORACLE_FAILED,
            diagnostics,
        )
        return OracleValidationResult(oracle_result, nothing, report, diagnostics)
    end

    if length(oracle_result.primal) != num_variables(problem)
        push!(
            diagnostics,
            "numerical oracle returned $(length(oracle_result.primal)) primal " *
            "entries; expected $(num_variables(problem))",
        )
        report = _failed_validation_report(NUMERICAL_ORACLE_FAILED, diagnostics)
        return OracleValidationResult(oracle_result, nothing, report, diagnostics)
    end

    hint = try
        NumericalPrimalHint(oracle_result.primal; objective = oracle_result.objective_value)
    catch error
        push!(diagnostics, "invalid numerical primal hint: $(sprint(showerror, error))")
        report = _failed_validation_report(NUMERICAL_ORACLE_FAILED, diagnostics)
        return OracleValidationResult(oracle_result, nothing, report, diagnostics)
    end

    recovery_options = if isnothing(options)
        RecoveryOptions(;
            max_denominator = max_denominator,
            atol = atol,
            rtol = rtol,
            objective_atol = objective_atol,
            objective_rtol = objective_rtol,
        )
    else
        options
    end

    certificate = try
        recover_primal_certificate(
            problem,
            hint,
            recovery_options;
            return_diagnostics = false,
        )::ExactPrimalCertificate
    catch error
        if error isa RationalRecoveryError
            recovery = error.diagnostics
            push!(
                diagnostics,
                "rational recovery failed at $(recovery.stage): $(recovery.message)",
            )
            report =
                isnothing(recovery.certificate_report) ?
                _failed_validation_report(recovery.status, diagnostics) :
                recovery.certificate_report
            return OracleValidationResult(oracle_result, nothing, report, diagnostics)
        end
        push!(
            diagnostics,
            "rational recovery threw an exception: $(sprint(showerror, error))",
        )
        report = _failed_validation_report(RECOVERY_FAILED, diagnostics)
        return OracleValidationResult(oracle_result, nothing, report, diagnostics)
    end

    report = check_certificate(problem, certificate; diagnostics = true)
    append!(diagnostics, report.diagnostics)
    report.ok || push!(diagnostics, "independent exact certificate checking failed")
    return OracleValidationResult(oracle_result, certificate, report, diagnostics)
end
