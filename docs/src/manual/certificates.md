# Certificates

`ExactPrimalCertificate` stores the rational point, package certificate version,
problem hash, exactification policy, objective value, and metadata.

`check_certificate(problem, certificate)` independently checks:

1. certificate version and problem hash;
2. exact affine equalities;
3. each cone block;
4. the exact objective value, when present.

Recovery and checking are deliberately separate. A checker never calls a numerical
solver or trusts recovery metadata.
