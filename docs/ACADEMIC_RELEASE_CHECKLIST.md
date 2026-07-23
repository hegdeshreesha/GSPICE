# Academic beta release checklist

Target: GSPICE 1.3.0-academic-beta

## Required before publishing

- [ ] Freeze scope: advertise only capabilities listed by `--capabilities`.
- [ ] Run clean Windows and Ubuntu CI from a fresh checkout.
- [ ] Run all CTests with SuiteSparse/KLU enabled and with the internal fallback.
- [ ] Run `DAE_AUDIT=1` decks across representative diode, BJT, and primitive
      MOS operating regions; archive the derivative and conservation reports.
- [ ] Run the state-preserving OSDI derivative/conservation audit across
      representative compact-model operating regions and supported ABIs.
- [ ] Measure convergence order for BE, trapezoidal, BDF2-5, and Adams2-5 on
      smooth analytic RLC/diode problems with fixed and changing step ratios.
- [ ] Compare OP/DC/TRAN/AC results on public RLC, diode, BJT, and OSDI decks
      against at least Xyce and ngspice; archive inputs, versions, tolerances,
      outputs, and comparison scripts.
- [ ] Add at least one independent numeric oracle for every advertised beta
      analysis; do not count a success-string regex as an accuracy oracle.
- [ ] Validate deterministic Monte Carlo repeatability for fixed seeds on both
      release platforms.
- [ ] Validate every bundled OSDI binary on both release platforms and record
      its compiler, ABI version, source/license, and SHA-256.
- [ ] Confirm `LICENSE`, `CITATION.cff`, dependency licenses, third-party model
      licenses, and PDK redistribution terms with the institution's release
      owner.
- [ ] Remove private/foundry/confidential decks, absolute paths, credentials,
      generated RAW/CSV files, and unlicensed model binaries from the archive.
- [ ] Publish `docs/LIMITATIONS.md` prominently in the release notes.
- [ ] State “academic beta; not for signoff” on the release page and binaries.
- [ ] Tag the exact tested commit; build artifacts from that tag, not a working
      directory.
- [ ] Generate checksums and malware-scan the archives.
- [ ] Have a second person reproduce the build and two reference comparisons.

## Recommended release contents

- Source archive and exact Git tag.
- Windows x64 and Ubuntu x86-64 binaries from CI.
- `LICENSE`, `NOTICE` if third-party notices become necessary, `CITATION.cff`,
  `README.md`, limitations, comparison, beta guide, and changelog.
- Public example decks plus expected outputs and tolerance-based validation.
- Dependency/model manifest with versions, origins, licenses, and hashes.

## Go/no-go criteria

Release only if all required items are complete, no supported deck can be made
to return a known false-success result, all CI tests pass, and the published
claims match `gspice --capabilities`. Otherwise delay the release or reduce the
advertised scope.

## After publishing

- Open a public issue tracker and security contact.
- Preserve release artifacts and reference results immutably.
- Triage correctness reports ahead of new features.
- Do not promote an experimental capability to beta until it has independent
  numeric oracles and cross-simulator validation.
