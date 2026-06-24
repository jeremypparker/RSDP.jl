using BenchmarkTools
using RSDP

const Q = Rational{BigInt}
const GRAM = let
    B = Q[1 2 0 1; 0 1 1 -1; 2 0 1 1]
    transpose(B) * B
end

@benchmark check_psd_exact($GRAM)
