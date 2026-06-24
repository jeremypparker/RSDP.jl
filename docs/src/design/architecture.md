# Architecture

The dependency direction is:

```text
exactification → exact affine algebra → cone checks
                         ↓                 ↓
                    recovery         certificate checker
                         ↘             ↙
                        exact certificate
```

MOI extraction and numerical oracles terminate at the exact problem or numerical-hint
boundary. They do not participate in checking. Optional ecosystem integrations must
not make solver or polynomial packages transitive dependencies of the checker.

The v0.1 scalar backend is `Rational{BigInt}` behind conversion and linear-algebra
functions. This is sufficient to validate the architecture before adding a Nemo,
sparse fraction-free, or modular backend.
