module RSDPHypatiaExt

using Hypatia
using RSDP
import MathOptInterface as MOI

const MOIU = MOI.Utilities

_float64(value::RSDP.ExactScalar) = Float64(value)

function _objective_constant(problem::RSDP.ExactConicProblem)
    value = get(problem.metadata, :objective_constant, zero(RSDP.ExactScalar))
    return _float64(RSDP.exactify(value; context = "objective constant"))
end

function _identity_function(
    variables::AbstractVector{MOI.VariableIndex},
    range::UnitRange{Int},
    coefficients::AbstractVector{Float64} = ones(Float64, length(range)),
)
    length(coefficients) == length(range) ||
        throw(DimensionMismatch("one coefficient is required per cone coordinate"))
    terms = MOI.VectorAffineTerm{Float64}[
        MOI.VectorAffineTerm(
            output,
            MOI.ScalarAffineTerm(coefficients[output], variables[column]),
        ) for (output, column) in enumerate(range)
    ]
    return MOI.VectorAffineFunction(terms, zeros(Float64, length(range)))
end

function _add_block!(
    model,
    variables,
    block::RSDP.ZeroConeBlock,
    first_column::Int,
    optimizer,
)
    range = first_column:(first_column+RSDP.dimension(block)-1)
    MOI.add_constraint(
        model,
        _identity_function(variables, range),
        MOI.Zeros(length(range)),
    )
    return last(range) + 1
end

function _add_block!(
    model,
    variables,
    block::RSDP.NonnegativeConeBlock,
    first_column::Int,
    optimizer,
)
    range = first_column:(first_column+RSDP.dimension(block)-1)
    MOI.add_constraint(
        model,
        _identity_function(variables, range),
        MOI.Nonnegatives(length(range)),
    )
    return last(range) + 1
end

function _add_block!(
    model,
    variables,
    block::RSDP.PSDTriangleConeBlock,
    first_column::Int,
    optimizer,
)
    range = first_column:(first_column+RSDP.dimension(block)-1)
    unscaled_set = MOI.PositiveSemidefiniteConeTriangle(block.side_dimension)
    if MOI.supports_constraint(
        optimizer,
        MOI.VectorAffineFunction{Float64},
        typeof(unscaled_set),
    )
        # Hypatia 0.9 supports the ordinary MOI triangle directly.
        function_value = _identity_function(variables, range)
        set = unscaled_set
    else
        # Hypatia 0.10 supports MOI's scaled PSD triangle. Scale the affine
        # function, not the RSDP variables, by sqrt(2) off diagonal. The
        # extracted primal therefore remains in RSDP's upper-triangle order:
        # (1,1), (1,2), (2,2), (1,3), ...
        scaling = collect(MOIU.SetDotScalingVector{Float64}(unscaled_set))
        function_value = _identity_function(variables, range, scaling)
        set = MOI.Scaled(unscaled_set)
    end
    MOI.add_constraint(model, function_value, set)
    return last(range) + 1
end

function _add_block!(
    model,
    variables,
    block::RSDP.ProductConeBlock,
    first_column::Int,
    optimizer,
)
    next_column = first_column
    for component in block.blocks
        next_column = _add_block!(model, variables, component, next_column, optimizer)
    end
    return next_column
end

function _unsupported_cone(block::RSDP.AbstractConeBlock)
    if block isa RSDP.ProductConeBlock
        for component in block.blocks
            unsupported = _unsupported_cone(component)
            isnothing(unsupported) || return unsupported
        end
        return nothing
    elseif block isa
           Union{RSDP.ZeroConeBlock, RSDP.NonnegativeConeBlock, RSDP.PSDTriangleConeBlock}
        return nothing
    end
    return block
end

function _instantiate(oracle::RSDP.HypatiaOracle)
    factory = oracle.optimizer_factory
    if isnothing(factory)
        return Hypatia.Optimizer()
    elseif factory isa MOI.OptimizerWithAttributes
        return MOI.instantiate(factory)
    elseif factory isa Type || factory isa Function
        return factory()
    end
    return factory
end

function _optimizer_attribute(key)
    if key isa MOI.AbstractOptimizerAttribute
        return key
    elseif key isa Symbol || key isa AbstractString
        return MOI.RawOptimizerAttribute(String(key))
    end
    throw(ArgumentError("unsupported optimizer attribute key $(repr(key))"))
end

function _safe_get(model, attribute, default)
    try
        return MOI.get(model, attribute)
    catch
        return default
    end
end

function RSDP.solve_oracle(
    problem::RSDP.ExactConicProblem,
    oracle::RSDP.HypatiaOracle;
    kwargs...,
)
    !isempty(kwargs) && return RSDP.NumericalOracleResult(
        RSDP.NUMERICAL_ORACLE_FAILED;
        raw_status = :unsupported_keyword,
        diagnostics = [
            "HypatiaOracle.solve_oracle does not accept keywords: " *
            join(string.(keys(kwargs)), ", "),
        ],
    )

    if !isnothing(problem.cones)
        unsupported = _unsupported_cone(problem.cones)
        if !isnothing(unsupported)
            return RSDP.NumericalOracleResult(
                RSDP.UNSUPPORTED_CONE;
                raw_status = :unsupported_cone,
                diagnostics = [
                    "HypatiaOracle does not support cone block " *
                    string(typeof(unsupported)),
                ],
            )
        end
    end

    optimizer = try
        _instantiate(oracle)
    catch error
        return RSDP.NumericalOracleResult(
            RSDP.NUMERICAL_ORACLE_FAILED;
            raw_status = :optimizer_instantiation_failed,
            raw_result = error,
            diagnostics = [
                "could not instantiate Hypatia optimizer: " * sprint(showerror, error),
            ],
        )
    end

    source = MOIU.Model{Float64}()
    variables = MOI.add_variables(source, RSDP.num_variables(problem))

    equality_terms = MOI.VectorAffineTerm{Float64}[]
    for row in axes(problem.A, 1), column in axes(problem.A, 2)
        coefficient = _float64(problem.A[row, column])
        iszero(coefficient) && continue
        push!(
            equality_terms,
            MOI.VectorAffineTerm(row, MOI.ScalarAffineTerm(coefficient, variables[column])),
        )
    end
    equality_function = MOI.VectorAffineFunction(equality_terms, -_float64.(problem.b))
    MOI.add_constraint(source, equality_function, MOI.Zeros(RSDP.num_constraints(problem)))

    if !isnothing(problem.cones)
        next_column = _add_block!(source, variables, problem.cones, 1, optimizer)
        next_column == RSDP.num_variables(problem) + 1 || return RSDP.NumericalOracleResult(
            RSDP.NUMERICAL_ORACLE_FAILED;
            raw_status = :internal_dimension_mismatch,
            diagnostics = ["cone translation did not consume every primal coordinate"],
        )
    end

    MOI.set(source, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    if !isnothing(problem.c)
        objective_terms = MOI.ScalarAffineTerm{Float64}[]
        for column in eachindex(problem.c)
            coefficient = _float64(problem.c[column])
            iszero(coefficient) && continue
            push!(objective_terms, MOI.ScalarAffineTerm(coefficient, variables[column]))
        end
        MOI.set(
            source,
            MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
            MOI.ScalarAffineFunction(objective_terms, _objective_constant(problem)),
        )
    end

    try
        MOI.set(optimizer, MOI.Silent(), oracle.silent)
        for (key, value) in oracle.attributes
            MOI.set(optimizer, _optimizer_attribute(key), value)
        end
    catch error
        return RSDP.NumericalOracleResult(
            RSDP.NUMERICAL_ORACLE_FAILED;
            raw_status = :optimizer_attribute_failed,
            raw_result = error,
            diagnostics = [
                "could not configure Hypatia optimizer: " * sprint(showerror, error),
            ],
        )
    end

    index_map = try
        MOI.copy_to(optimizer, source)
    catch error
        return RSDP.NumericalOracleResult(
            RSDP.NUMERICAL_ORACLE_FAILED;
            raw_status = :model_copy_failed,
            raw_result = error,
            diagnostics = [
                "could not copy the floating conic model to Hypatia: " *
                sprint(showerror, error),
            ],
        )
    end

    try
        MOI.optimize!(optimizer)
    catch error
        return RSDP.NumericalOracleResult(
            RSDP.NUMERICAL_ORACLE_FAILED;
            raw_status = :optimizer_exception,
            raw_result = optimizer,
            diagnostics = ["Hypatia failed: $(sprint(showerror, error))"],
        )
    end

    termination = _safe_get(optimizer, MOI.TerminationStatus(), MOI.OTHER_ERROR)
    primal_status = _safe_get(optimizer, MOI.PrimalStatus(), MOI.NO_SOLUTION)
    result_count = _safe_get(optimizer, MOI.ResultCount(), 0)
    raw_status = (
        termination = termination,
        primal_status = primal_status,
        result_count = result_count,
        solver_status = _safe_get(optimizer, MOI.RawStatusString(), ""),
    )
    diagnostics = String[
        "Hypatia termination status: $termination",
        "Hypatia primal status: $primal_status",
    ]

    usable_status =
        primal_status == MOI.FEASIBLE_POINT || primal_status == MOI.NEARLY_FEASIBLE_POINT
    if result_count < 1 || !usable_status
        push!(diagnostics, "Hypatia returned no usable primal point")
        return RSDP.NumericalOracleResult(
            RSDP.NUMERICAL_ORACLE_FAILED;
            raw_status = raw_status,
            raw_result = optimizer,
            diagnostics = diagnostics,
        )
    end

    mapped_variables = MOI.VariableIndex[index_map[variable] for variable in variables]
    primal = try
        Float64[
            MOI.get(optimizer, MOI.VariablePrimal(), variable) for
            variable in mapped_variables
        ]
    catch error
        push!(
            diagnostics,
            "could not extract Hypatia primal values: $(sprint(showerror, error))",
        )
        return RSDP.NumericalOracleResult(
            RSDP.NUMERICAL_ORACLE_FAILED;
            raw_status = raw_status,
            raw_result = optimizer,
            diagnostics = diagnostics,
        )
    end
    if !all(isfinite, primal)
        push!(diagnostics, "Hypatia primal vector contains non-finite values")
        return RSDP.NumericalOracleResult(
            RSDP.NUMERICAL_ORACLE_FAILED;
            raw_status = raw_status,
            raw_result = optimizer,
            diagnostics = diagnostics,
        )
    end

    objective_value = if isnothing(problem.c)
        nothing
    else
        value = _safe_get(optimizer, MOI.ObjectiveValue(), nothing)
        value isa Real && isfinite(value) ? Float64(value) : nothing
    end
    return RSDP.NumericalOracleResult(
        RSDP.NUMERICAL_SOLVED_NOT_VALIDATED;
        primal = primal,
        objective_value = objective_value,
        raw_status = raw_status,
        raw_result = optimizer,
        diagnostics = diagnostics,
    )
end

end
