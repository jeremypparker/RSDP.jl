using Test
import MathOptInterface as MOI

const MOIU = MOI.Utilities

function _scalar_affine(terms, constant = 0)
    Q = Rational{BigInt}
    return MOI.ScalarAffineFunction{Q}(
        MOI.ScalarAffineTerm{Q}[
            MOI.ScalarAffineTerm(Q(coefficient), variable) for
            (coefficient, variable) in terms
        ],
        Q(constant),
    )
end

function _vector_affine(terms, constants)
    Q = Rational{BigInt}
    return MOI.VectorAffineFunction{Q}(
        MOI.VectorAffineTerm{Q}[
            MOI.VectorAffineTerm(output, MOI.ScalarAffineTerm(Q(coefficient), variable)) for
            (output, coefficient, variable) in terms
        ],
        Q.(constants),
    )
end

@testset "strict MathOptInterface extraction" begin
    Q = Rational{BigInt}
    model = MOIU.Model{Q}()
    x = MOI.add_variables(model, 2)

    # Add constraints in deliberately noncanonical order. Extraction order must
    # depend only on supported family and MOI index, not dictionary iteration.
    psd_variables = MOI.add_constraint(
        model,
        MOI.VectorOfVariables([x[1], x[2], x[1]]),
        MOI.PositiveSemidefiniteConeTriangle(2),
    )
    nonnegative_variables =
        MOI.add_constraint(model, MOI.VectorOfVariables([x[2], x[1]]), MOI.Nonnegatives(2))
    equality = MOI.add_constraint(
        model,
        _scalar_affine([(1, x[1]), (2, x[2])], 1),
        MOI.EqualTo(Q(4)),
    )
    psd_affine = MOI.add_constraint(
        model,
        _vector_affine([(1, 1, x[1]), (2, 1, x[2]), (3, 1, x[1])], [0, 0, 1]),
        MOI.PositiveSemidefiniteConeTriangle(2),
    )
    zero_variables = MOI.add_constraint(model, MOI.VectorOfVariables([x[2]]), MOI.Zeros(1))
    nonnegative_affine = MOI.add_constraint(
        model,
        _vector_affine([(1, 1, x[1]), (2, 2, x[2])], [1, -1]),
        MOI.Nonnegatives(2),
    )
    zero_affine = MOI.add_constraint(
        model,
        _vector_affine([(1, 1, x[1]), (1, -1, x[2]), (2, 3, x[2])], [2, -1]),
        MOI.Zeros(2),
    )
    objective = _scalar_affine([(2, x[1]), (-3, x[2])], 5)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Q}}(), objective)

    extracted = RSDP.extract_moi(model)

    @test extracted.variables == x
    @test extracted.positive_columns == 1:2
    @test extracted.negative_columns == 3:4
    @test size(extracted.A) == (14, 14)
    @test extracted.b[1:4] == Q[3, -2, 1, 0]
    @test extracted.A[1:4, 1:4] == Q[
        1 2 -1 -2
        1 -1 -1 1
        0 3 0 -3
        0 1 0 -1
    ]
    @test extracted.c == Q[2, -3, -2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    @test extracted.objective_constant == Q(5)
    @test extracted.objective_sense == MOI.MIN_SENSE
    @test RSDP.equivalent(extracted) === extracted.problem
    @test RSDP.num_constraints(extracted) == 14
    @test RSDP.num_variables(extracted) == 14

    @test extracted.block_columns == [1:2, 3:4, 5:6, 7:8, 9:11, 12:14]
    @test extracted.block_sources == [
        :moi_variable_positive_parts,
        :moi_variable_negative_parts,
        nonnegative_affine,
        nonnegative_variables,
        psd_affine,
        psd_variables,
    ]
    @test extracted.blocks == [
        RSDP.NonnegativeConeBlock(2),
        RSDP.NonnegativeConeBlock(2),
        RSDP.NonnegativeConeBlock(2),
        RSDP.NonnegativeConeBlock(2),
        RSDP.PSDTriangleConeBlock(2),
        RSDP.PSDTriangleConeBlock(2),
    ]
    @test extracted.row_sources == [
        (equality, 1),
        (zero_affine, 1),
        (zero_affine, 2),
        (zero_variables, 1),
        (nonnegative_affine, 1),
        (nonnegative_affine, 2),
        (nonnegative_variables, 1),
        (nonnegative_variables, 2),
        (psd_affine, 1),
        (psd_affine, 2),
        (psd_affine, 3),
        (psd_variables, 1),
        (psd_variables, 2),
        (psd_variables, 3),
    ]

    # Affine conic coordinates satisfy y = f(u).
    @test extracted.A[5:6, 1:6] == Q[
        -1 0 1 0 1 0
        0 -2 0 2 0 1
    ]
    @test extracted.b[5:6] == Q[1, -1]
    @test extracted.A[9:11, 1:11] == Q[
        -1 0 1 0 0 0 0 0 1 0 0
        0 -1 0 1 0 0 0 0 0 1 0
        -1 0 1 0 0 0 0 0 0 0 1
    ]
    @test extracted.b[9:11] == Q[0, 0, 1]

    coordinates = Q[4, 3, 1, 2, 4, 5, 3, 3, 3, 1, 4, 3, 1, 3]
    @test RSDP.recover_moi_variables(extracted, coordinates) == Q[3, 1]
    @test_throws DimensionMismatch RSDP.recover_moi_variables(extracted, Q[1])
end

@testset "objective normalization and exactness policy" begin
    Q = Rational{BigInt}
    exact_model = MOIU.Model{Q}()
    x = MOI.add_variable(exact_model)
    MOI.set(exact_model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.set(
        exact_model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Q}}(),
        _scalar_affine([(2, x)], 3),
    )
    maximum = RSDP.MOIExtractedProblem(exact_model)
    @test maximum.c == Q[-2, 2]
    @test maximum.objective_constant == Q(-3)
    @test maximum.objective_sense == MOI.MAX_SENSE
    certificate = RSDP.make_primal_certificate(maximum, Q[2, 0])
    report = RSDP.check_certificate(maximum, certificate)
    @test report.ok
    @test report.computed_objective == -7

    float_model = MOIU.Model{Float64}()
    y = MOI.add_variable(float_model)
    MOI.add_constraint(
        float_model,
        MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(0.1, y)], 0.2),
        MOI.EqualTo(0.3),
    )
    MOI.set(float_model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(
        float_model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(0.5, y)], 0.0),
    )
    @test_throws RSDP.InexactDataError RSDP.extract_moi(float_model)

    decimal = RSDP.extract_moi(float_model; policy = RSDP.DecimalStringInexact())
    @test decimal.A == Q[1//10 -1//10]
    @test decimal.b == Q[1//10]
    @test decimal.c == Q[1//2, -1//2]

    binary = RSDP.extract_moi(float_model; policy = RSDP.RationalizeInexact(0.0))
    @test binary.A[1, 1] == Q(BigInt(3602879701896397), BigInt(36028797018963968))
end

@testset "strictly reject unsupported MOI forms" begin
    Q = Rational{BigInt}
    constraint_model = MOIU.Model{Q}()
    x = MOI.add_variable(constraint_model)
    MOI.add_constraint(constraint_model, _scalar_affine([(1, x)]), MOI.LessThan(Q(1)))
    @test_throws RSDP.InvalidProblemError RSDP.extract_moi(constraint_model)

    objective_model = MOIU.Model{Q}()
    y = MOI.add_variable(objective_model)
    MOI.set(objective_model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(objective_model, MOI.ObjectiveFunction{MOI.VariableIndex}(), y)
    @test_throws RSDP.InvalidProblemError RSDP.extract_moi(objective_model)
end
