const MOIU = MOI.Utilities

"""
    ValidationStatusAttribute()

MOI optimizer attribute returning RSDP's exact validation status.
"""
struct ValidationStatusAttribute <: MOI.AbstractOptimizerAttribute end

"""
    CertificateAttribute()

MOI optimizer attribute returning the recovered exact primal certificate, or
`nothing` when no certificate was recovered.
"""
struct CertificateAttribute <: MOI.AbstractOptimizerAttribute end

"""
    ValidationReportAttribute()

MOI optimizer attribute returning the independent exact validation report, or
`nothing` before validation reaches the certificate checker.
"""
struct ValidationReportAttribute <: MOI.AbstractOptimizerAttribute end

"""
    OracleResultAttribute()

MOI optimizer attribute returning the numerical-oracle result, or `nothing`
when extraction failed before the oracle was called.
"""
struct OracleResultAttribute <: MOI.AbstractOptimizerAttribute end

"""
    DiagnosticsAttribute()

MOI optimizer attribute returning validation and extraction diagnostics.
"""
struct DiagnosticsAttribute <: MOI.AbstractOptimizerAttribute end

"""
    MaxDenominatorAttribute()

Settable MOI optimizer attribute controlling rational recovery's denominator
bound.
"""
struct MaxDenominatorAttribute <: MOI.AbstractOptimizerAttribute end

"""
    NumericalOracleAttribute()

Settable MOI optimizer attribute containing the numerical hint generator.
"""
struct NumericalOracleAttribute <: MOI.AbstractOptimizerAttribute end

"""
    ExactificationPolicyAttribute()

Settable MOI optimizer attribute controlling conversion of model data to exact
rationals. The default is [`ErrorOnInexact`](@ref).
"""
struct ExactificationPolicyAttribute <: MOI.AbstractOptimizerAttribute end

for Attribute in (
    ValidationStatusAttribute,
    CertificateAttribute,
    ValidationReportAttribute,
    OracleResultAttribute,
    DiagnosticsAttribute,
)
    @eval MOI.is_set_by_optimize(::$Attribute) = true
    @eval MOIU.map_indices(::Any, ::$Attribute, value) = value
end

"""
    Optimizer(; oracle=HypatiaOracle(), max_denominator=10^8,
                exactification_policy=ErrorOnInexact(), silent=true)

A thin MathOptInterface optimizer that validates primal feasibility through
RSDP's existing MOI extraction, numerical-oracle, rational-recovery, and exact
certificate-checking pipeline.

The optimizer does not prove optimality. Hypatia remains optional and is used
only through the numerical-oracle extension.
"""
mutable struct Optimizer <: MOI.AbstractOptimizer
    model::MOIU.Model{ExactScalar}
    oracle::AbstractNumericalOracle
    max_denominator::BigInt
    exactification_policy::AbstractInexactPolicy
    validation_result::Union{Nothing, OracleValidationResult}
    extraction_result::Union{Nothing, MOIExtractedProblem}
    variable_values::Dict{MOI.VariableIndex, ExactScalar}
    status::ValidationStatus
    raw_status::String
    silent::Bool
    diagnostics_log::Vector{String}
end

function Optimizer(;
    oracle::AbstractNumericalOracle = HypatiaOracle(),
    max_denominator::Integer = big(10)^8,
    exactification_policy::AbstractInexactPolicy = ErrorOnInexact(),
    silent::Bool = true,
)
    max_denominator > 0 || throw(ArgumentError("max_denominator must be positive"))
    return Optimizer(
        MOIU.Model{ExactScalar}(),
        oracle,
        BigInt(max_denominator),
        exactification_policy,
        nothing,
        nothing,
        Dict{MOI.VariableIndex, ExactScalar}(),
        NOT_SOLVED,
        "RSDP has not been run",
        silent,
        String[],
    )
end

function _clear_result!(optimizer::Optimizer)
    optimizer.validation_result = nothing
    optimizer.extraction_result = nothing
    empty!(optimizer.variable_values)
    optimizer.status = NOT_SOLVED
    optimizer.raw_status = "RSDP has not been run"
    empty!(optimizer.diagnostics_log)
    return nothing
end

MOI.supports_incremental_interface(::Optimizer) = true
MOI.copy_to(destination::Optimizer, source::MOI.ModelLike) =
    MOIU.default_copy_to(destination, source)

function MOI.empty!(optimizer::Optimizer)
    MOI.empty!(optimizer.model)
    _clear_result!(optimizer)
    return nothing
end

MOI.is_empty(optimizer::Optimizer) = MOI.is_empty(optimizer.model)
MOI.add_variable(optimizer::Optimizer) = MOI.add_variable(optimizer.model)
MOI.add_variables(optimizer::Optimizer, count::Int) =
    MOI.add_variables(optimizer.model, count)

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.ScalarAffineFunction{ExactScalar}},
    ::Type{MOI.EqualTo{ExactScalar}},
)
    return true
end

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{F},
    ::Type{S},
) where {
    F <: Union{MOI.VectorAffineFunction{ExactScalar}, MOI.VectorOfVariables},
    S <: Union{
        MOI.Zeros,
        MOI.Nonnegatives,
        MOI.PositiveSemidefiniteConeTriangle,
    },
}
    return true
end

MOI.supports_constraint(
    ::Optimizer,
    ::Type{F},
    ::Type{S},
) where {F <: MOI.AbstractFunction, S <: MOI.AbstractSet} = false

function MOI.add_constraint(
    optimizer::Optimizer,
    func::F,
    set::S,
) where {F <: MOI.AbstractFunction, S <: MOI.AbstractSet}
    MOI.supports_constraint(optimizer, F, S) ||
        throw(MOI.UnsupportedConstraint{F, S}())
    return MOI.add_constraint(optimizer.model, func, set)
end

MOI.is_valid(optimizer::Optimizer, index) = MOI.is_valid(optimizer.model, index)
MOI.delete(optimizer::Optimizer, index) = MOI.delete(optimizer.model, index)
MOI.modify(optimizer::Optimizer, index, change) =
    MOI.modify(optimizer.model, index, change)

function MOI.supports(optimizer::Optimizer, attribute::MOI.AbstractModelAttribute)
    if attribute isa MOI.ObjectiveFunction
        return attribute isa MOI.ObjectiveFunction{MOI.ScalarAffineFunction{ExactScalar}}
    end
    return MOI.supports(optimizer.model, attribute)
end

MOI.supports(optimizer::Optimizer, attribute::MOI.AbstractVariableAttribute, index_type) =
    MOI.supports(optimizer.model, attribute, index_type)
MOI.supports(optimizer::Optimizer, attribute::MOI.AbstractConstraintAttribute, index_type) =
    MOI.supports(optimizer.model, attribute, index_type)

MOI.get(optimizer::Optimizer, attribute::MOI.AbstractModelAttribute) =
    MOI.get(optimizer.model, attribute)
MOI.get(
    optimizer::Optimizer,
    attribute::MOI.AbstractVariableAttribute,
    variable::MOI.VariableIndex,
) = MOI.get(optimizer.model, attribute, variable)
MOI.get(
    optimizer::Optimizer,
    attribute::MOI.AbstractConstraintAttribute,
    constraint::MOI.ConstraintIndex,
) = MOI.get(optimizer.model, attribute, constraint)

MOI.set(
    optimizer::Optimizer,
    attribute::MOI.AbstractModelAttribute,
    value,
) = MOI.set(optimizer.model, attribute, value)
MOI.set(
    optimizer::Optimizer,
    attribute::MOI.AbstractVariableAttribute,
    variable::MOI.VariableIndex,
    value,
) = MOI.set(optimizer.model, attribute, variable, value)
MOI.set(
    optimizer::Optimizer,
    attribute::MOI.AbstractConstraintAttribute,
    constraint::MOI.ConstraintIndex,
    value,
) = MOI.set(optimizer.model, attribute, constraint, value)

MOI.get(optimizer::Optimizer, ::Type{MOI.VariableIndex}, name::String) =
    MOI.get(optimizer.model, MOI.VariableIndex, name)
MOI.get(optimizer::Optimizer, ::Type{MOI.ConstraintIndex}, name::String) =
    MOI.get(optimizer.model, MOI.ConstraintIndex, name)
MOI.get(
    optimizer::Optimizer,
    ::Type{MOI.ConstraintIndex{F, S}},
    name::String,
) where {F, S} = MOI.get(optimizer.model, MOI.ConstraintIndex{F, S}, name)

MOI.supports(::Optimizer, ::MOI.Silent) = true
MOI.get(optimizer::Optimizer, ::MOI.Silent) = optimizer.silent
MOI.set(optimizer::Optimizer, ::MOI.Silent, value::Bool) = (optimizer.silent = value)
MOI.get(::Optimizer, ::MOI.SolverName) = "RSDP validated primal feasibility"

for Attribute in (
    ValidationStatusAttribute,
    CertificateAttribute,
    ValidationReportAttribute,
    OracleResultAttribute,
    DiagnosticsAttribute,
    MaxDenominatorAttribute,
    NumericalOracleAttribute,
    ExactificationPolicyAttribute,
)
    @eval MOI.supports(::Optimizer, ::$Attribute) = true
end

MOI.get(optimizer::Optimizer, ::ValidationStatusAttribute) = optimizer.status
MOI.get(optimizer::Optimizer, ::CertificateAttribute) =
    isnothing(optimizer.validation_result) ? nothing :
    optimizer.validation_result.certificate
MOI.get(optimizer::Optimizer, ::ValidationReportAttribute) =
    isnothing(optimizer.validation_result) ? nothing :
    optimizer.validation_result.report
MOI.get(optimizer::Optimizer, ::OracleResultAttribute) =
    isnothing(optimizer.validation_result) ? nothing :
    optimizer.validation_result.oracle_result
MOI.get(optimizer::Optimizer, ::DiagnosticsAttribute) =
    copy(optimizer.diagnostics_log)
MOI.get(optimizer::Optimizer, ::MaxDenominatorAttribute) =
    optimizer.max_denominator
MOI.get(optimizer::Optimizer, ::NumericalOracleAttribute) = optimizer.oracle
MOI.get(optimizer::Optimizer, ::ExactificationPolicyAttribute) =
    optimizer.exactification_policy

function MOI.set(
    optimizer::Optimizer,
    ::MaxDenominatorAttribute,
    value::Integer,
)
    value > 0 || throw(ArgumentError("max_denominator must be positive"))
    optimizer.max_denominator = BigInt(value)
    return nothing
end

function MOI.set(
    optimizer::Optimizer,
    ::NumericalOracleAttribute,
    oracle::AbstractNumericalOracle,
)
    optimizer.oracle = oracle
    return nothing
end

function MOI.set(
    optimizer::Optimizer,
    ::ExactificationPolicyAttribute,
    policy::AbstractInexactPolicy,
)
    optimizer.exactification_policy = policy
    return nothing
end

function _extraction_failure_status(error)
    if error isa InexactDataError || error isa InexactConversionError
        return EXACTIFICATION_REQUIRED
    elseif error isa InvalidProblemError
        return UNSUPPORTED_MODEL
    end
    return UNSUPPORTED_MODEL
end

function _configured_oracle(optimizer::Optimizer)
    oracle = optimizer.oracle
    if oracle isa HypatiaOracle
        return HypatiaOracle(
            optimizer_factory = oracle.optimizer_factory,
            attributes = oracle.attributes,
            silent = optimizer.silent,
        )
    end
    return oracle
end

function MOI.optimize!(optimizer::Optimizer)
    _clear_result!(optimizer)
    extracted = try
        extract_moi(
            optimizer.model;
            policy = optimizer.exactification_policy,
        )
    catch error
        optimizer.status = _extraction_failure_status(error)
        message = "MOI extraction failed: $(sprint(showerror, error))"
        push!(optimizer.diagnostics_log, message)
        optimizer.raw_status = message
        return nothing
    end
    optimizer.extraction_result = extracted

    result = validate_with_oracle(
        extracted.problem,
        _configured_oracle(optimizer);
        max_denominator = optimizer.max_denominator,
    )
    optimizer.validation_result = result
    append!(optimizer.diagnostics_log, result.diagnostics)
    optimizer.status = result.report.ok ? VALIDATED_PRIMAL_FEASIBLE :
                       result.report.status
    optimizer.raw_status = result.report.ok ?
                           "RSDP validated primal feasibility" :
                           "RSDP validation failed: " *
                           (isempty(result.diagnostics) ? string(optimizer.status) :
                            join(result.diagnostics, "; "))

    if !isnothing(result.certificate)
        values = recover_moi_variables(extracted, result.certificate.x)
        for (variable, value) in zip(extracted.variables, values)
            optimizer.variable_values[variable] = value
        end
    end
    return nothing
end

function MOI.get(optimizer::Optimizer, ::MOI.TerminationStatus)
    if optimizer.status == VALIDATED_PRIMAL_FEASIBLE
        return MOI.LOCALLY_SOLVED
    elseif optimizer.status in (EXACTIFICATION_REQUIRED, UNSUPPORTED_MODEL, UNSUPPORTED_CONE)
        return MOI.INVALID_MODEL
    elseif optimizer.status == NOT_SOLVED
        return MOI.OPTIMIZE_NOT_CALLED
    elseif optimizer.status == NUMERICAL_ORACLE_FAILED
        return MOI.NUMERICAL_ERROR
    end
    return MOI.OTHER_ERROR
end

MOI.get(optimizer::Optimizer, ::MOI.PrimalStatus) =
    optimizer.status == VALIDATED_PRIMAL_FEASIBLE ? MOI.FEASIBLE_POINT : MOI.NO_SOLUTION
MOI.get(::Optimizer, ::MOI.DualStatus) = MOI.NO_SOLUTION
MOI.get(optimizer::Optimizer, ::MOI.ResultCount) =
    optimizer.status == VALIDATED_PRIMAL_FEASIBLE ? 1 : 0
MOI.get(optimizer::Optimizer, ::MOI.RawStatusString) = optimizer.raw_status

function MOI.get(
    optimizer::Optimizer,
    attribute::MOI.VariablePrimal,
    variable::MOI.VariableIndex,
)
    MOI.check_result_index_bounds(optimizer, attribute)
    return optimizer.variable_values[variable]
end

function MOI.get(optimizer::Optimizer, attribute::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(optimizer, attribute)
    sense = MOI.get(optimizer.model, MOI.ObjectiveSense())
    sense == MOI.FEASIBILITY_SENSE && return zero(ExactScalar)
    objective_type = MOI.get(optimizer.model, MOI.ObjectiveFunctionType())
    objective = MOI.get(
        optimizer.model,
        MOI.ObjectiveFunction{objective_type}(),
    )
    value = exactify(objective.constant)
    for term in objective.terms
        value += exactify(term.coefficient) * optimizer.variable_values[term.variable]
    end
    return value
end

"""
    validation_status(model)

Return RSDP's validation status from an optimizer or JuMP/MOI model.
"""
validation_status(model) = MOI.get(model, ValidationStatusAttribute())

"""
    certificate(model)

Return the exact primal certificate from an optimizer or JuMP/MOI model.
"""
certificate(model) = MOI.get(model, CertificateAttribute())

"""
    validation_report(model)

Return the exact certificate-checking report from an optimizer or JuMP/MOI
model.
"""
validation_report(model) = MOI.get(model, ValidationReportAttribute())

"""
    oracle_result(model)

Return the numerical-oracle result from an optimizer or JuMP/MOI model.
"""
oracle_result(model) = MOI.get(model, OracleResultAttribute())

"""
    diagnostics(model)

Return validation diagnostics from an optimizer or JuMP/MOI model.
"""
diagnostics(model) = MOI.get(model, DiagnosticsAttribute())
