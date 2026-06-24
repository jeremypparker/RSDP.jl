# Design spikes

## Spike 1: hand-built conic data

A deterministic test constructs a rational product-cone problem, perturbs a feasible
point by hand, recovers affine coordinates, reconstructs an exact point, and invokes
the independent checker.

## Spike 2: MOI extraction

The same problem is represented in an MOI model. Extraction is compared with the
hand-built problem under deterministic ordering.

## Spike 3: SOS and WeightedSOSCone

The supported early path is bridged MOI validation. SOS-aware reconstruction remains
conditional on stable public metadata; version-dependent access belongs in an
extension and must have version guards.

## Spike 4: boundary face

A rank-deficient PSD example exercises boundary diagnostics. Automatic facial
reduction is deferred until exact exposing-vector and face-basis certificates can be
replayed by the checker.
