# Literature context

RSDP's affine-space rounding follows the exact-validation pattern used by
Dostert, de Laat, and Moustrou for exact SDP bounds. Its SOS recovery direction is
informed by Peyrl--Parrilo and by Kaltofen--Li--Yang--Zhi. Laplagne's work motivates
certified facial reduction for degenerate SOS problems.

Scheiderer's counterexamples are a central limitation: a real SOS or feasible real SDP
need not admit a rational certificate. Davis--Papp-style rational dual WSOS
certificates are future work, as are rational-function certificates and formal-proof
exports.

The implemented v0.1 claim is deliberately smaller: exact checking of rational primal
conic certificates, plus heuristic recovery whose output must pass that checker.
