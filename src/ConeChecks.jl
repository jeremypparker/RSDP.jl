export ConeMembershipDiagnostic,
    ConeMembershipResult, check_cone_membership, check_cone, check_cones, in_cone

"""
    ConeMembershipDiagnostic

A diagnostic emitted by [`check_cone_membership`](@ref). `path` identifies a
component within nested product cones, and `index` identifies a packed
coordinate when applicable.
"""
struct ConeMembershipDiagnostic
    code::Symbol
    message::String
    path::Vector{Int}
    index::Union{Nothing, Int}
    psd_result::Union{Nothing, PSDCheckResult}
end

"""
    ConeMembershipResult

Structured exact cone-membership result.
"""
struct ConeMembershipResult
    is_member::Bool
    diagnostics::Vector{ConeMembershipDiagnostic}
end

function Base.show(io::IO, result::ConeMembershipResult)
    print(
        io,
        "ConeMembershipResult(",
        result.is_member ? "member" : "not a member",
        ", diagnostics=",
        length(result.diagnostics),
        ")",
    )
end

function _cone_diagnostic(
    code::Symbol,
    message::String,
    path::Vector{Int};
    index::Union{Nothing, Int} = nothing,
    psd_result::Union{Nothing, PSDCheckResult} = nothing,
)
    return ConeMembershipDiagnostic(code, message, copy(path), index, psd_result)
end

function _cone_exact_values(values::AbstractVector, path::Vector{Int})
    exact = Vector{ExactRational}(undef, length(values))
    for index in eachindex(values)
        try
            exact[index] = _exact_rational(values[index])
        catch error
            error isa ArgumentError || rethrow()
            diagnostic = _cone_diagnostic(
                :inexact_input,
                "coordinate $index has non-exact type $(typeof(values[index])); " *
                "use integers or rationals",
                path;
                index = Int(index),
            )
            return nothing, diagnostic
        end
    end
    return exact, nothing
end

function _dimension_diagnostic(
    block::AbstractConeBlock,
    values::AbstractVector,
    path::Vector{Int},
)
    expected = dimension(block)
    length(values) == expected && return nothing
    return _cone_diagnostic(
        :dimension_mismatch,
        "cone block requires $expected coordinates, got $(length(values))",
        path,
    )
end

function _check_cone_membership(
    block::ZeroConeBlock,
    values::AbstractVector,
    path::Vector{Int},
)
    mismatch = _dimension_diagnostic(block, values, path)
    mismatch === nothing || return ConeMembershipResult(false, [mismatch])
    exact, diagnostic = _cone_exact_values(values, path)
    diagnostic === nothing || return ConeMembershipResult(false, [diagnostic])

    bad_index = findfirst(value -> !iszero(value), exact)
    if bad_index === nothing
        return ConeMembershipResult(true, ConeMembershipDiagnostic[])
    end
    failure = _cone_diagnostic(
        :nonzero_in_zero_cone,
        "zero-cone coordinate $bad_index is nonzero",
        path;
        index = bad_index,
    )
    return ConeMembershipResult(false, [failure])
end

function _check_cone_membership(
    block::NonnegativeConeBlock,
    values::AbstractVector,
    path::Vector{Int},
)
    mismatch = _dimension_diagnostic(block, values, path)
    mismatch === nothing || return ConeMembershipResult(false, [mismatch])
    exact, diagnostic = _cone_exact_values(values, path)
    diagnostic === nothing || return ConeMembershipResult(false, [diagnostic])

    bad_index = findfirst(value -> value < 0, exact)
    if bad_index === nothing
        return ConeMembershipResult(true, ConeMembershipDiagnostic[])
    end
    failure = _cone_diagnostic(
        :negative_coordinate,
        "nonnegative-cone coordinate $bad_index is negative",
        path;
        index = bad_index,
    )
    return ConeMembershipResult(false, [failure])
end

function _check_cone_membership(
    block::PSDTriangleConeBlock,
    values::AbstractVector,
    path::Vector{Int},
)
    mismatch = _dimension_diagnostic(block, values, path)
    mismatch === nothing || return ConeMembershipResult(false, [mismatch])
    exact, diagnostic = _cone_exact_values(values, path)
    diagnostic === nothing || return ConeMembershipResult(false, [diagnostic])

    matrix = triangle_to_matrix(exact, block)
    psd_result = check_psd(matrix)
    if psd_result.is_psd
        return ConeMembershipResult(true, ConeMembershipDiagnostic[])
    end
    failure = _cone_diagnostic(
        :not_positive_semidefinite,
        psd_result.diagnostic.message,
        path;
        index = psd_result.diagnostic.index,
        psd_result = psd_result,
    )
    return ConeMembershipResult(false, [failure])
end

function _check_cone_membership(
    block::ProductConeBlock,
    values::AbstractVector,
    path::Vector{Int},
)
    mismatch = _dimension_diagnostic(block, values, path)
    mismatch === nothing || return ConeMembershipResult(false, [mismatch])

    diagnostics = ConeMembershipDiagnostic[]
    offset = 0
    for (component_index, component) in enumerate(block.blocks)
        component_dimension = dimension(component)
        component_values = view(values, (offset+1):(offset+component_dimension))
        component_path = vcat(path, component_index)
        result = _check_cone_membership(component, component_values, component_path)
        append!(diagnostics, result.diagnostics)
        offset += component_dimension
    end
    return ConeMembershipResult(isempty(diagnostics), diagnostics)
end

"""
    check_cone_membership(block, values) -> ConeMembershipResult

Check membership in a cone block using exact integer/rational arithmetic.
Product-cone diagnostics carry a path of one-based component indices.
"""
function check_cone_membership(block::AbstractConeBlock, values::AbstractVector)
    return _check_cone_membership(block, values, Int[])
end

"""
Alias for [`check_cone_membership`](@ref).
"""
check_cone(block::AbstractConeBlock, values::AbstractVector) =
    check_cone_membership(block, values)

"""
    check_cones(block_or_blocks, values)

Check one cone block, or a tuple/vector of blocks interpreted as a product.
"""
check_cones(block::AbstractConeBlock, values::AbstractVector) =
    check_cone_membership(block, values)
check_cones(blocks::Tuple, values::AbstractVector) =
    check_cone_membership(ProductConeBlock(blocks), values)
check_cones(blocks::AbstractVector{<:AbstractConeBlock}, values::AbstractVector) =
    check_cone_membership(ProductConeBlock(blocks), values)

"""
Return only the Boolean outcome of [`check_cone_membership`](@ref).
"""
in_cone(block::AbstractConeBlock, values::AbstractVector) =
    check_cone_membership(block, values).is_member

"""
    check_cones(problem, values)

Check a problem's full product-cone membership exactly.
"""
function check_cones(problem::ExactConicProblem, values::AbstractVector)
    isnothing(problem.cones) && throw(InvalidProblemError("problem has no cone metadata"))
    return check_cone_membership(problem.cones, values)
end
