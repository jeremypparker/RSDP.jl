# MOI extraction

The extractor accepts a deliberately small MathOptInterface subset: scalar and vector
affine equalities, `Zeros`, `Nonnegatives`,
`PositiveSemidefiniteConeTriangle`, and scalar affine objectives.

Variable and block order are deterministic. Unsupported function/set pairs fail
before a partial certificate problem is returned. Inexact coefficients are rejected
unless the caller supplies an exactification policy.

JuMP and SumOfSquares models should be bridged to this supported conic subset before
extraction. The bridge result, not an unbridged symbolic model, is what the conic
certificate validates.
