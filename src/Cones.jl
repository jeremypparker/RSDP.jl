"""
    ExactRational

The scalar type used by RSDP's exact cone and positive-semidefinite checks.
Floating-point values are deliberately not converted to this type implicitly.
"""
const ExactRational = Rational{BigInt}

export ExactRational,
    AbstractConeBlock,
    ZeroConeBlock,
    NonnegativeConeBlock,
    PSDTriangleConeBlock,
    ProductConeBlock,
    triangle_dimension,
    triangle_side_dimension,
    dimension,
    ambient_dimension,
    matrix_dimension,
    triangle_to_matrix,
    matrix_to_triangle,
    reconstruct_triangle

"""
    AbstractConeBlock

Supertype for finite-dimensional cone blocks.
"""
abstract type AbstractConeBlock end

function _checked_dimension(value::Integer, name::AbstractString)
    value >= 0 || throw(ArgumentError("$name must be nonnegative, got $value"))
    try
        return Int(value)
    catch error
        error isa InexactError || rethrow()
        throw(ArgumentError("$name does not fit in Int: $value"))
    end
end

"""
    ZeroConeBlock(dimension)

The zero cone `{0}` in the indicated ambient dimension.
"""
struct ZeroConeBlock <: AbstractConeBlock
    dimension::Int

    function ZeroConeBlock(dimension::Integer)
        return new(_checked_dimension(dimension, "cone dimension"))
    end
end

"""
    NonnegativeConeBlock(dimension)

The nonnegative orthant in the indicated ambient dimension.
"""
struct NonnegativeConeBlock <: AbstractConeBlock
    dimension::Int

    function NonnegativeConeBlock(dimension::Integer)
        return new(_checked_dimension(dimension, "cone dimension"))
    end
end

"""
    PSDTriangleConeBlock(side_dimension)

The cone of `side_dimension × side_dimension` positive-semidefinite matrices,
represented by their upper triangles in column-major order:
`(1,1), (1,2), (2,2), (1,3), ...`.
"""
struct PSDTriangleConeBlock <: AbstractConeBlock
    side_dimension::Int

    function PSDTriangleConeBlock(side_dimension::Integer)
        return new(_checked_dimension(side_dimension, "PSD side dimension"))
    end
end

"""
    ProductConeBlock(blocks...)
    ProductConeBlock(blocks)

The Cartesian product of zero or more cone blocks.
"""
struct ProductConeBlock{B<:Tuple} <: AbstractConeBlock
    blocks::B

    function ProductConeBlock(blocks::B) where {B<:Tuple}
        all(block -> block isa AbstractConeBlock, blocks) ||
            throw(ArgumentError("every product component must be an AbstractConeBlock"))
        return new{B}(blocks)
    end
end

ProductConeBlock(blocks::AbstractConeBlock...) = ProductConeBlock(blocks)
ProductConeBlock(blocks::AbstractVector{<:AbstractConeBlock}) =
    ProductConeBlock(tuple(blocks...))

"""
    triangle_dimension(side_dimension)

Return the number of entries in one triangle of a square matrix.
"""
function triangle_dimension(side_dimension::Integer)
    n = _checked_dimension(side_dimension, "triangle side dimension")
    try
        return Base.checked_mul(n, Base.checked_add(n, 1)) ÷ 2
    catch error
        error isa OverflowError || rethrow()
        throw(ArgumentError("triangle dimension overflows Int for side dimension $n"))
    end
end

"""
    triangle_side_dimension(packed_dimension)

Recover the side dimension from a triangular packed-vector length. Throw an
`ArgumentError` when the length is not triangular.
"""
function triangle_side_dimension(packed_dimension::Integer)
    m = _checked_dimension(packed_dimension, "packed triangle dimension")
    discriminant = BigInt(1) + BigInt(8) * BigInt(m)
    root = isqrt(discriminant)
    root * root == discriminant || throw(ArgumentError("$m is not a triangular number"))
    isodd(root) || throw(ArgumentError("$m is not a triangular number"))
    n = (root - 1) ÷ 2
    triangle_dimension(n) == m || throw(ArgumentError("$m is not a triangular number"))
    return Int(n)
end

"""
    dimension(block)

Return the packed ambient dimension of a cone block.
"""
dimension(block::Union{ZeroConeBlock,NonnegativeConeBlock}) = block.dimension
dimension(block::PSDTriangleConeBlock) = triangle_dimension(block.side_dimension)

function dimension(block::ProductConeBlock)
    total = 0
    for component in block.blocks
        try
            total = Base.checked_add(total, dimension(component))
        catch error
            error isa OverflowError || rethrow()
            throw(ArgumentError("product cone dimension overflows Int"))
        end
    end
    return total
end

"""Alias for [`dimension`](@ref)."""
ambient_dimension(block::AbstractConeBlock) = dimension(block)

"""Return the matrix side dimension of a PSD triangle block."""
matrix_dimension(block::PSDTriangleConeBlock) = block.side_dimension

Base.length(block::AbstractConeBlock) = dimension(block)

_exact_rational(value::ExactRational) = value
_exact_rational(value::Integer) = BigInt(value) // BigInt(1)
_exact_rational(value::Rational{T}) where {T<:Integer} =
    BigInt(numerator(value)) // BigInt(denominator(value))

function _exact_rational(value)
    throw(ArgumentError(
        "expected an integer or rational value for an exact check, got $(typeof(value))",
    ))
end

"""
    triangle_to_matrix(values[, side_dimension])

Reconstruct an exact symmetric matrix from an upper-triangular packed vector.
The vector order is column-major within the upper triangle:
`(1,1), (1,2), (2,2), (1,3), ...`. Integer and rational entries are converted
to `Rational{BigInt}`; inexact entries are rejected.
"""
function triangle_to_matrix(values::AbstractVector, side_dimension::Integer)
    n = _checked_dimension(side_dimension, "triangle side dimension")
    expected = triangle_dimension(n)
    length(values) == expected || throw(DimensionMismatch(
        "packed triangle has length $(length(values)); expected $expected " *
        "for side dimension $n",
    ))

    matrix = Matrix{ExactRational}(undef, n, n)
    packed_index = 1
    for column in 1:n
        for row in 1:column
            value = _exact_rational(values[packed_index])
            matrix[row, column] = value
            matrix[column, row] = value
            packed_index += 1
        end
    end
    return matrix
end

triangle_to_matrix(values::AbstractVector) =
    triangle_to_matrix(values, triangle_side_dimension(length(values)))
triangle_to_matrix(values::AbstractVector, block::PSDTriangleConeBlock) =
    triangle_to_matrix(values, block.side_dimension)

"""
    matrix_to_triangle(matrix)

Pack a symmetric exact matrix using the convention of
[`triangle_to_matrix`](@ref). Nonsquare, nonsymmetric, and inexact matrices are
rejected.
"""
function matrix_to_triangle(matrix::AbstractMatrix)
    rows, columns = size(matrix)
    rows == columns || throw(DimensionMismatch(
        "cannot pack a nonsquare $(rows) × $(columns) matrix",
    ))

    result = Vector{ExactRational}(undef, triangle_dimension(rows))
    packed_index = 1
    for column in 1:columns
        for row in 1:column
            upper = _exact_rational(matrix[row, column])
            lower = _exact_rational(matrix[column, row])
            upper == lower || throw(ArgumentError(
                "matrix is not symmetric at indices ($row, $column) and ($column, $row)",
            ))
            result[packed_index] = upper
            packed_index += 1
        end
    end
    return result
end

"""Alias for [`triangle_to_matrix`](@ref)."""
reconstruct_triangle(args...) = triangle_to_matrix(args...)
