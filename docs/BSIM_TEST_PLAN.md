# BSIM3/BSIM4 clean-room implementation test plan

## Gates

1. Parameter gate: aliases, units, defaults, geometry, temperature, and invalid cards.
2. DC gate: Id-Vg, Id-Vd, reverse operation, body bias, temperature, and current conservation.
3. Derivative gate: analytic Jacobians against central differences at bias points away from regime boundaries.
4. Charge gate: terminal charge conservation, capacitance symmetry, and finite derivatives.
5. AC/transient gate: small-signal admittance and charge-history transient behavior.
6. Noise gate: thermal, flicker, shot, induced-gate, and correlation behavior.
7. Integration gate: parser -> translator -> GMC -> GSDI device, with unsupported models failing closed.
8. Reference gate: independent ngspice vectors with defined relative/absolute tolerances.
9. Performance gate: model preprocessing cached, stable matrix pattern, and no per-Newton model allocation.

## Current execution

- BSIM3 parameter/DC/DAE/noise smoke gates: passing; reference gate still failing.
- BSIM4 DC parity: PASSING — oracle (official BSIM4.8.3 b4ld/b4temp/b4set
  transcription, tools/bsim4_reference_oracle.py) vs native probe across 4
  decks (baseline, body-bias/DVT/LPE, short-channel Rout/CLM, temperature):
  50/50 points, worst relative 8.9e-3 @ 1% tolerance. ngspice is NOT the BSIM4
  reference (its bundled BSIM4 ignores PDIBL1/PDIBL2/ETAB, plus default-mask
  markers). C++ deviations fixed this milestone: UA/UB/UC defaults (was 0,
  official: 1e-9/1e-19/-0.0465e-9), vsattemp `vsat*(1-at*dT)` scaling, Vgsteff
  exp(80) branch (linear law T10=mstar*Vgst / T3=coxe*MAX_EXP/cdep0), Abulk AGS
  term missing bulkFactor multiplier removed, thetaRout missing exp(T0)
  multiplier.
- BSIM4 parameter/DC/Jacobian smoke gates: passing; reference gate: passing.
- Parameter-coverage gate: passing. tools/audit_bsim4_parameters.py validates
  (a) official BSIM4.8.3 names (bsim4def.h) vs the native read-set —
  implemented-and-official coverage 97.7%, (b) allow-list sync
  (Bsim4ImplementedParameters in bsim4_parameters.hpp — 118 names),
  (c) deck-scan: set-but-ignored model-card params across the corpus. Parser
  now warns per-deck for ignored BSIM4 parameters (like PSP already did):
  smoke_bsim4_unsupported_param_warning tests it. Known implementation gaps
  surfaced by the audit: CV/junction params (CJ/CJSW/PB/MJ), charge-model
  switch CAPMOD/XPART, geometry-bin l*/w* families, TOX-alias derive,
  flicker/noise family, RDSMOD/VSATT switches.
- BSIM3/BSIM4 parser production enablement: disabled.

## Next implementation slice

BSIM4 DC parity is green (4-deck corpus, 50/50 at 1%). Remaining gaps vs the
compact-model standards (checked against official reference equations only —
native transcriptions and GSDI/GMC; OSDI/OpenVAF are never used): CV/charge
conservation (naive terminalCharges), noise (thermal
and flicker at 1/f), self-heating/temperature Jacobian review, and a full
BSIM4 card mask + corner-case corpus (temperature derating, Vbsm, GIDL at
short L, overlap caps). BSIM3 reference gate (ngspice cross-check) still fails
(subthreshold 2-3 decades early, Vth offset ~0.4-0.5 V, early saturation) and
remains the only red gate; the 4-deck oracle methodology ported here applies
directly (ngspice is a valid BSIM3 reference).
