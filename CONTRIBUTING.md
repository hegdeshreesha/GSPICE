# Contributing

Thank you for improving GSPICE. Correctness and honest capability reporting take
priority over feature count.

## Clean-room implementation

GSPICE may learn requirements and numerical behavior from public simulators,
standards, and academic literature, but contributions must be independently
implemented. Do not copy source code, pseudocode, comments, tests, internal
names, or file organization from VACASK, Xyce, ngspice, or other projects whose
license is incompatible with GSPICE. See
[`docs/DAE_ARCHITECTURE.md`](docs/DAE_ARCHITECTURE.md) for the design and review
rules used by the charge-model migration.

1. Create a focused branch and keep unrelated changes separate.
2. Add a regression test for every behavior change. Numerical features should
   include tolerance-based reference values, ideally from an independently
   implemented analytic result or a documented Xyce/ngspice comparison.
3. Never silently ignore an active device, source, directive, or failed solve.
4. Mark unvalidated analyses/models experimental in code, documentation, and
   `--capabilities`.
5. Build and run `ctest --test-dir build -C Release --output-on-failure`.
6. Update `docs/LIMITATIONS.md` when a limitation is added, narrowed, or removed.

By contributing, you agree that your contribution is licensed under
Apache-2.0 and that you have the right to submit it. Do not contribute
confidential PDK content, restricted foundry models, or code copied from
incompatibly licensed simulators.
