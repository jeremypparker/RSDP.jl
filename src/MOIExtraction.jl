import MathOptInterface as MOI

export MOIExtractedProblem,
    extract_moi,
    extract_moi_problem,
    equivalent,
    recover_moi_variables

"""
    MOIExtractedProblem

An exact, flat cone-coordinate representation extracted from a
MathOptInterface model.

`problem` stores the normalized minimization problem `A * x = b` with
objective `c'x + objective_constant`. `cone` partitions every coordinate of
`x`. Original MOI variables are unrestricted, so variable `u[j]` is represented
exactly as

```math
u[j] = x[positive_columns[j]] - x[negative_columns[j]],
```

where both coordinate groups are nonnegative cone blocks. Each affine
nonnegative or PSD constraint gets a fresh cone block and equality rows linking
that block to its affine function. `Zeros` and scalar `EqualTo` constraints
become equality rows directly.

Variables are ordered by increasing `MOI.VariableIndex.value`. Constraint
families have the canonical order documented by [`extract_moi`](@ref), and
indices within each family are ordered by increasing constraint-index value.
"""
struct MOIExtractedProblem
    problem::ExactConicProblem
    cone::ProductConeBlock
    objective_constant::ExactScalar
    objective_sense::MOI.OptimizationSense
    variables::Vector{MOI.VariableIndex}
    positive_columns::UnitRange{Int}
    negative_columns::UnitRange{Int}
    block_columns::Vector{UnitRange{Int}}
    block_sources::Vector{Any}
    row_sources::Vector{Tuple{Any,Int}}

    function MOIExtractedProblem(
        problem::ExactConicProblem,
        cone::ProductConeBlock,
        objective_constant::ExactScalar,
        objective_sense::MOI.OptimizationSense,
        variables::Vector{MOI.VariableIndex},
        positive_columns::UnitRange{Int},
        negative_columns::UnitRange{Int},
        block_columns::Vector{UnitRange{Int}},
        block_sources::Vector{Any},
        row_sources::Vector{Tuple{Any,Int}},
    )
        dimension(cone) == num_variables(problem) ||
            throw(
                InvalidProblemError(
                    "cone dimension $(dimension(cone)) does not match " *
                    "$(num_variables(problem)) flat coordinates",
                ),
            )
        length(block_columns) == length(cone.blocks) ||
            throw(InvalidProblemError("one column range is required per cone block"))
        length(block_sources) == length(cone.blocks) ||
            throw(InvalidProblemError("one source is required per cone block"))
        length(row_sources) == num_constraints(problem) ||
            throw(InvalidProblemError("one source is required per equality row"))
        return new(
            problem,
            cone,
            objective_constant,
            objective_sense,
            copy(variables),
            positive_columns,
            negative_columns,
            copy(block_columns),
            copy(block_sources),
            copy(row_sources),
        )
    end
end

function Base.getproperty(extracted::MOIExtractedProblem, name::Symbol)
    if name === :A || name === :b || name === :c
        return getproperty(getfield(extracted, :problem), name)
    elseif name === :cone_blocks || name === :blocks
        return collect(getfield(extracted, :cone).blocks)
    end
    return getfield(extracted, name)
end

function Base.propertynames(::MOIExtractedProblem, private::Bool=false)
    public = (
        :problem,
        :A,
        :b,
        :c,
        :cone,
        :cone_blocks,
        :blocks,
        :objective_constant,
        :objective_sense,
        :variables,
        :positive_columns,
        :negative_columns,
        :block_columns,
        :block_sources,
        :row_sources,
    )
    return private ? public : public
end

"""
    equivalent(extracted::MOIExtractedProblem)

Return the core [`ExactConicProblem`](@ref) represented by `extracted`.
Cone-block and MOI-origin metadata remain available on `extracted`.
"""
equivalent(extracted::MOIExtractedProblem) = extracted.problem

num_constraints(extracted::MOIExtractedProblem) =
    num_constraints(extracted.problem)
num_variables(extracted::MOIExtractedProblem) = num_variables(extracted.problem)

"""
    recover_moi_variables(extracted, x)

Recover original MOI-variable values from a flat cone-coordinate vector `x`.
This only applies the exact free-variable split; it does not check `A*x == b`
or cone membership.
"""
function recover_moi_variables(
    extracted::MOIExtractedProblem,
    x::AbstractVector,
)
    length(x) == num_variables(extracted.problem) ||
        throw(
            DimensionMismatch(
                "expected $(num_variables(extracted.problem)) cone coordinates, " *
                "got $(length(x))",
            ),
        )
    return x[extracted.positive_columns] - x[extracted.negative_columns]
end

struct _MOIConstraintRecord
    rank::Int
    index::Any
    func::Any
    set::Any
end

function _constraint_rank(F::Type, S::Type)
    if F <: MOI.ScalarAffineFunction && S <: MOI.EqualTo
        return 1
    elseif F <: MOI.VectorAffineFunction && S <: MOI.Zeros
        return 2
    elseif F === MOI.VectorOfVariables && S <: MOI.Zeros
        return 3
    elseif F <: MOI.VectorAffineFunction && S <: MOI.Nonnegatives
        return 4
    elseif F === MOI.VectorOfVariables && S <: MOI.Nonnegatives
        return 5
    elseif F <: MOI.VectorAffineFunction &&
           S <: MOI.PositiveSemidefiniteConeTriangle
        return 6
    elseif F === MOI.VectorOfVariables &&
           S <: MOI.PositiveSemidefiniteConeTriangle
        return 7
    end
    return 0
end

_constraint_index_value(index) = getfield(index, :value)

function _collect_constraints(model::MOI.ModelLike)
    records = _MOIConstraintRecord[]
    types_present = MOI.get(model, MOI.ListOfConstraintTypesPresent())
    ordered_types = sort(
        collect(types_present);
        by=pair -> begin
            F, S = pair
            rank = _constraint_rank(F, S)
            return (rank == 0 ? typemax(Int) : rank, string(F), string(S))
        end,
    )
    for (F, S) in ordered_types
        rank = _constraint_rank(F, S)
        rank != 0 ||
            throw(
                InvalidProblemError(
                    "unsupported MOI constraint type $F-in-$S; supported forms are " *
                    "ScalarAffineFunction-in-EqualTo and VectorAffineFunction or " *
                    "VectorOfVariables in Zeros, Nonnegatives, or " *
                    "PositiveSemidefiniteConeTriangle",
                ),
            )
        indices = sort(
            collect(MOI.get(model, MOI.ListOfConstraintIndices{F,S}()));
            by=_constraint_index_value,
        )
        for index in indices
            func = MOI.get(model, MOI.ConstraintFunction(), index)
            set = MOI.get(model, MOI.ConstraintSet(), index)
            _validate_constraint_dimension(index, func, set)
            push!(records, _MOIConstraintRecord(rank, index, func, set))
        end
    end
    return records
end

_function_dimension(::MOI.ScalarAffineFunction) = 1
_function_dimension(func::MOI.VectorAffineFunction) = length(func.constants)
_function_dimension(func::MOI.VectorOfVariables) = length(func.variables)

function _validate_constraint_dimension(index, func, set)
    function_dimension = _function_dimension(func)
    set_dimension = set isa MOI.EqualTo ? 1 : MOI.dimension(set)
    function_dimension == set_dimension ||
        throw(
            InvalidProblemError(
                "constraint $index has function dimension $function_dimension " *
                "and set dimension $set_dimension",
            ),
        )
    if func isa MOI.VectorAffineFunction
        for term in func.terms
            1 <= term.output_index <= function_dimension ||
                throw(
                    InvalidProblemError(
                        "constraint $index has vector-affine output index " *
                        "$(term.output_index) outside 1:$function_dimension",
                    ),
                )
        end
    end
    return nothing
end

_record_dimension(record::_MOIConstraintRecord) = _function_dimension(record.func)
_is_conic_record(record::_MOIConstraintRecord) = record.rank >= 4

function _exact_moi(value, policy::AbstractInexactPolicy, context)
    return exactify(value, policy; context=context)
end

function _variable_position(
    positions::Dict{MOI.VariableIndex,Int},
    variable::MOI.VariableIndex,
    context,
)
    position = get(positions, variable, 0)
    position != 0 ||
        throw(
            InvalidProblemError(
                "$context references variable $variable, which is not in " *
                "MOI.ListOfVariableIndices",
            ),
        )
    return position
end

function _add_original_coefficient!(
    A::Matrix{ExactScalar},
    row::Int,
    variable::MOI.VariableIndex,
    coefficient::ExactScalar,
    positions::Dict{MOI.VariableIndex,Int},
    variable_count::Int,
    context,
)
    position = _variable_position(positions, variable, context)
    A[row, position] += coefficient
    A[row, variable_count + position] -= coefficient
    return nothing
end

function _add_scalar_equality!(
    A,
    b,
    row,
    record,
    positions,
    variable_count,
    policy,
)
    func = record.func
    context = "constraint $(record.index)"
    for term in func.terms
        coefficient = _exact_moi(
            term.coefficient,
            policy,
            "$context coefficient for $(term.variable)",
        )
        _add_original_coefficient!(
            A,
            row,
            term.variable,
            coefficient,
            positions,
            variable_count,
            context,
        )
    end
    rhs = _exact_moi(record.set.value, policy, "$context EqualTo value")
    constant = _exact_moi(func.constant, policy, "$context constant")
    b[row] = rhs - constant
    return nothing
end

function _add_vector_zero_rows!(
    A,
    b,
    first_row,
    record,
    positions,
    variable_count,
    policy,
)
    context = "constraint $(record.index)"
    func = record.func
    if func isa MOI.VectorAffineFunction
        for term in func.terms
            row = first_row + term.output_index - 1
            scalar_term = term.scalar_term
            coefficient = _exact_moi(
                scalar_term.coefficient,
                policy,
                "$context output $(term.output_index) coefficient for " *
                "$(scalar_term.variable)",
            )
            _add_original_coefficient!(
                A,
                row,
                scalar_term.variable,
                coefficient,
                positions,
                variable_count,
                context,
            )
        end
        for output in eachindex(func.constants)
            b[first_row + output - 1] =
                -_exact_moi(func.constants[output], policy, "$context constant[$output]")
        end
    else
        for (output, variable) in enumerate(func.variables)
            _add_original_coefficient!(
                A,
                first_row + output - 1,
                variable,
                one(ExactScalar),
                positions,
                variable_count,
                context,
            )
        end
    end
    return nothing
end

function _add_conic_link_rows!(
    A,
    b,
    first_row,
    first_column,
    record,
    positions,
    variable_count,
    policy,
)
    context = "constraint $(record.index)"
    func = record.func
    dimension = _record_dimension(record)
    for output in 1:dimension
        A[first_row + output - 1, first_column + output - 1] = one(ExactScalar)
    end
    if func isa MOI.VectorAffineFunction
        for term in func.terms
            row = first_row + term.output_index - 1
            scalar_term = term.scalar_term
            coefficient = -_exact_moi(
                scalar_term.coefficient,
                policy,
                "$context output $(term.output_index) coefficient for " *
                "$(scalar_term.variable)",
            )
            _add_original_coefficient!(
                A,
                row,
                scalar_term.variable,
                coefficient,
                positions,
                variable_count,
                context,
            )
        end
        for output in eachindex(func.constants)
            b[first_row + output - 1] =
                _exact_moi(func.constants[output], policy, "$context constant[$output]")
        end
    else
        for (output, variable) in enumerate(func.variables)
            _add_original_coefficient!(
                A,
                first_row + output - 1,
                variable,
                -one(ExactScalar),
                positions,
                variable_count,
                context,
            )
        end
    end
    return nothing
end

function _cone_block(set::MOI.Nonnegatives)
    return NonnegativeConeBlock(MOI.dimension(set))
end

function _cone_block(set::MOI.PositiveSemidefiniteConeTriangle)
    return PSDTriangleConeBlock(set.side_dimension)
end

function _extract_objective!(
    c,
    model,
    positions,
    variable_count,
    policy,
)
    sense = MOI.get(model, MOI.ObjectiveSense())
    if sense == MOI.FEASIBILITY_SENSE
        return zero(ExactScalar), sense
    end
    objective_type = MOI.get(model, MOI.ObjectiveFunctionType())
    objective_type <: MOI.ScalarAffineFunction ||
        throw(
            InvalidProblemError(
                "unsupported MOI objective type $objective_type; expected a " *
                "ScalarAffineFunction",
            ),
        )
    objective = MOI.get(model, MOI.ObjectiveFunction{objective_type}())
    multiplier = if sense == MOI.MIN_SENSE
        one(ExactScalar)
    elseif sense == MOI.MAX_SENSE
        -one(ExactScalar)
    else
        throw(InvalidProblemError("unsupported MOI objective sense $sense"))
    end
    for term in objective.terms
        coefficient =
            multiplier * _exact_moi(
                term.coefficient,
                policy,
                "objective coefficient for $(term.variable)",
            )
        position =
            _variable_position(positions, term.variable, "objective")
        c[position] += coefficient
        c[variable_count + position] -= coefficient
    end
    constant =
        multiplier * _exact_moi(objective.constant, policy, "objective constant")
    return constant, sense
end

"""
    extract_moi(model; policy=ErrorOnInexact())

Extract the supported part of an `MOI.ModelLike` into exact cone-coordinate
standard form. The supported constraint families, in deterministic output
order, are:

1. scalar affine in `EqualTo`;
2. vector affine, then `VectorOfVariables`, in `Zeros`;
3. vector affine, then `VectorOfVariables`, in `Nonnegatives`;
4. vector affine, then `VectorOfVariables`, in
   `PositiveSemidefiniteConeTriangle`.

Original variables are ordered by increasing MOI index and split into positive
and negative parts. Nonnegative and PSD constraints receive fresh cone
coordinates in the order above. This avoids free cone coordinates while
remaining exactly equivalent to the MOI model.

All coefficients, constants, and equality right-hand sides pass through
[`exactify`](@ref). Therefore floating-point model data is rejected by default.
Pass an explicit policy such as `RationalizeInexact()` or
`DecimalStringInexact()` to accept it.

Scalar affine objectives are normalized to minimization. Maximization
objectives are negated, including `objective_constant`; `objective_sense`
records the original sense.
"""
function extract_moi(
    model::MOI.ModelLike;
    policy::AbstractInexactPolicy=DEFAULT_INEXACT_POLICY,
)
    variables = sort(
        collect(MOI.get(model, MOI.ListOfVariableIndices()));
        by=variable -> variable.value,
    )
    length(unique(variables)) == length(variables) ||
        throw(InvalidProblemError("MOI.ListOfVariableIndices contains duplicates"))
    variable_count = length(variables)
    positions = Dict(variable => index for (index, variable) in enumerate(variables))
    records = _collect_constraints(model)

    equality_rows = sum(
        _record_dimension(record) for record in records if !_is_conic_record(record);
        init=0,
    )
    conic_rows = sum(
        _record_dimension(record) for record in records if _is_conic_record(record);
        init=0,
    )
    conic_columns = conic_rows
    row_count = equality_rows + conic_rows
    column_count = 2 * variable_count + conic_columns

    A = zeros(ExactScalar, row_count, column_count)
    b = zeros(ExactScalar, row_count)
    c = zeros(ExactScalar, column_count)

    positive_columns = 1:variable_count
    negative_columns = (variable_count + 1):(2 * variable_count)
    blocks = AbstractConeBlock[]
    block_columns = UnitRange{Int}[]
    block_sources = Any[]
    if variable_count > 0
        push!(blocks, NonnegativeConeBlock(variable_count))
        push!(block_columns, positive_columns)
        push!(block_sources, :moi_variable_positive_parts)
        push!(blocks, NonnegativeConeBlock(variable_count))
        push!(block_columns, negative_columns)
        push!(block_sources, :moi_variable_negative_parts)
    end

    row_sources = Tuple{Any,Int}[]
    row = 1
    column = 2 * variable_count + 1
    for record in records
        record_dimension = _record_dimension(record)
        if record.rank == 1
            _add_scalar_equality!(
                A,
                b,
                row,
                record,
                positions,
                variable_count,
                policy,
            )
        elseif record.rank <= 3
            _add_vector_zero_rows!(
                A,
                b,
                row,
                record,
                positions,
                variable_count,
                policy,
            )
        else
            range = column:(column + record_dimension - 1)
            push!(blocks, _cone_block(record.set))
            push!(block_columns, range)
            push!(block_sources, record.index)
            _add_conic_link_rows!(
                A,
                b,
                row,
                column,
                record,
                positions,
                variable_count,
                policy,
            )
            column += record_dimension
        end
        for output in 1:record_dimension
            push!(row_sources, (record.index, output))
        end
        row += record_dimension
    end
    row == row_count + 1 ||
        throw(InvalidProblemError("internal MOI row-count mismatch"))
    column == column_count + 1 ||
        throw(InvalidProblemError("internal MOI column-count mismatch"))

    objective_constant, objective_sense =
        _extract_objective!(c, model, positions, variable_count, policy)
    cone = ProductConeBlock(blocks)
    problem = ExactConicProblem(
        A,
        b,
        cone;
        c=c,
        policy=policy,
        metadata=Dict(
            :source => :MathOptInterface,
            :objective_constant => objective_constant,
            :objective_sense => objective_sense,
        ),
    )
    return MOIExtractedProblem(
        problem,
        cone,
        objective_constant,
        objective_sense,
        variables,
        positive_columns,
        negative_columns,
        block_columns,
        block_sources,
        row_sources,
    )
end

extract_moi(model::MOI.ModelLike, policy::AbstractInexactPolicy) =
    extract_moi(model; policy=policy)

"""Compatibility alias for [`extract_moi`](@ref)."""
extract_moi_problem(args...; kwargs...) = extract_moi(args...; kwargs...)

"""Construct an extracted problem directly from an MOI model."""
MOIExtractedProblem(model::MOI.ModelLike; kwargs...) = extract_moi(model; kwargs...)
