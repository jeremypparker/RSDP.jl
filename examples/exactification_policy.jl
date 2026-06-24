using RSDP

@assert exactify(0.1, DecimalStringInexact()) == 1 // 10

try
    exactify(0.1)
catch error
    @assert error isa InexactDataError
end
