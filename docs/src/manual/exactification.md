# Exactification

The default `ErrorOnInexact()` policy rejects binary floating-point model data.
Integers and rational values convert exactly.

`RationalizeInexact` is an explicit approximation policy with a tolerance and
denominator bound. `DecimalStringInexact` interprets the printed decimal value as an
exact decimal rational. These policies answer different mathematical questions and
are recorded in certificate metadata.

Never use exactification to conceal uncertainty in source data. If coefficients were
measured or rounded, the resulting rational problem is a chosen surrogate problem.
