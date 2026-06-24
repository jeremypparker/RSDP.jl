"""
    SolveStatus

Status of an exact feasibility computation.
"""
@enum SolveStatus::UInt8 begin
    UNKNOWN = 0
    FEASIBLE = 1
    INFEASIBLE = 2
end

"""
    ExactAffineStatus

Alias for [`SolveStatus`](@ref), retained to make affine-result signatures
self-documenting.
"""
const ExactAffineStatus = SolveStatus

"""
Return `true` when `status` is `FEASIBLE`.
"""
is_feasible(status::SolveStatus) = status == FEASIBLE

"""
Return `true` when `status` is `INFEASIBLE`.
"""
is_infeasible(status::SolveStatus) = status == INFEASIBLE

"""
Return `true` when `status` is `UNKNOWN`.
"""
is_unknown(status::SolveStatus) = status == UNKNOWN

"""
    ValidationStatus

Status of numerical hint generation, exact recovery, or certificate validation.
Only statuses beginning with `VALIDATED_` represent exact mathematical claims.
"""
@enum ValidationStatus::UInt8 begin
    NOT_SOLVED = 0
    NUMERICAL_SOLVED_NOT_VALIDATED
    VALIDATED_PRIMAL_FEASIBLE
    VALIDATED_DUAL_FEASIBLE
    VALIDATED_DUAL_BOUND
    VALIDATED_OPTIMAL
    VALIDATED_INFEASIBLE
    EXACTIFICATION_REQUIRED
    UNSUPPORTED_CONE
    UNSUPPORTED_OBJECTIVE
    UNSUPPORTED_MODEL
    INCONSISTENT_AFFINE_SYSTEM
    RECOVERY_FAILED
    RECOVERY_FAILED_AFFINE
    RECOVERY_FAILED_CONE
    RECOVERY_FAILED_DENOMINATOR_LIMIT
    RECOVERY_FAILED_WITH_BOUNDARY_DIAGNOSTIC
    FACIAL_REDUCTION_REQUIRED
    CERTIFIED_FACE_REDUCTION_FAILED
    NON_RATIONAL_FACE_SUSPECTED
    NUMERICAL_ORACLE_FAILED
    CERTIFICATE_CHECK_FAILED
end
