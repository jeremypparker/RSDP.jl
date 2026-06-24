using BenchmarkTools
using RSDP

const Q = Rational{BigInt}
const A_AFFINE = Q[
    1 1 0 0
    0 1 1 0
    0 0 1 1
]
const B_AFFINE = Q[1, 1, 1]

@benchmark exact_affine_space($A_AFFINE, $B_AFFINE)
