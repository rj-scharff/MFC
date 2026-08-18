# MFC Research-Fork Agent Instructions

## Fork purpose

This repository is the MFC host-code fork for the sibling research project:

`../MFCNucleationCavitation`

The canonical governing strategy is:

`../MFCNucleationCavitation/docs/strategy.md`

Read that document before performing substantive work related to nucleation,
phase change, Euler–Euler bubbles, six-equation coupling, or energy accounting.

If the sibling strategy cannot be accessed, stop and report that limitation
rather than reconstructing the strategy from memory.

This file adds research-fork constraints. Continue to respect MFC's own
contribution, formatting, testing, and source-organization conventions.

## Current branch and phase

The authorized research branch is:

`feature/nucleation-cavitation-audit`

The current project phase is:

> Phase 0 — source audit and baseline reproduction.

Do not perform research work on `master`.

Do not treat the name of this feature branch as permission to implement the
eventual nucleation model.

## Default operating mode

The default mode is read-only source audit.

Permitted activities include:

- reading and searching source files;
- tracing configuration checks and model-selection logic;
- tracing state allocation and variable indexing;
- tracing source-term and relaxation call ordering;
- identifying existing examples and regression tests;
- building MFC;
- running existing examples and tests;
- recording exact commands and outputs in the companion repository;
- adding minimal read-only diagnostics only when the user explicitly requests a
  code change.

## Changes prohibited during Phase 0

Do not perform any of the following without explicit user authorization:

- implement nucleation;
- add a source to bubble number density or population moments;
- change the PDE state vector or `num_pdes`;
- alter conserved-to-primitive or primitive-to-conserved semantics;
- modify pressure-relaxation physics;
- modify temperature or chemical-potential relaxation physics;
- change bubble-wall dynamics;
- change existing bubble mass-transfer behavior;
- add or replace an equation of state;
- change WENO, TENO, HLL, HLLC, limiters, or reconstruction defaults;
- alter axisymmetric or cylindrical geometry algorithms;
- update regression golden files;
- weaken an input validation check merely to make a combined configuration run;
- perform broad cleanup, renaming, formatting, or refactoring;
- update dependencies or toolchains;
- open or submit an upstream pull request.

If an incompatibility is discovered, document it before proposing a fix.

## Rules for explicitly authorized instrumentation

When the user explicitly authorizes an audit instrumentation change:

1. Make the smallest possible change.
2. Keep all existing behavior unchanged when the diagnostic is disabled.
3. Prefer a new, narrowly scoped file or routine over invasive edits.
4. Put new behavior behind a default-off input or compile-time control when
   appropriate.
5. Do not alter existing numerical results merely to simplify instrumentation.
6. Do not change golden files to make a new result pass.
7. Run the narrowest relevant build and regression tests.
8. Record the exact commands and outputs.
9. Report compiler warnings and failures; do not hide them.
10. Verify `git diff --check` before presenting the change.

## Required audit questions

Source inspection shall establish:

1. Whether the six-equation carrier and Euler–Euler bubble paths are mutually
   permitted by preprocessing and input validation.
2. Whether phase change and Euler–Euler bubbles can be active simultaneously.
3. Whether Strang splitting supports the selected carrier or only a single
   resolved liquid phase.
4. The ownership and storage of:
   - liquid and vapor mass;
   - volume fraction;
   - phasic internal energy;
   - mixture energy;
   - bubble pressure;
   - bubble vapor mass;
   - bubble number density;
   - class variables;
   - moment variables.
5. The exact order of:
   - finite-volume transport;
   - Runge–Kutta stages;
   - primitive conversion;
   - pressure relaxation;
   - thermal and chemical relaxation;
   - bubble source integration;
   - limiting;
   - reinitialization or projection.
6. Whether bubble micro-inertial, thermal, surface, or vapor energy is already
   represented in any existing energy variable.
7. Whether the moment path can accept a birth source without violating the
   inversion assumptions.
8. The exact requirements of MFC's axisymmetric mode.

## Evidence discipline

For every conclusion, cite the exact:

- source path;
- module;
- subroutine or function;
- variable or configuration parameter;
- relevant branch or conditional;
- pinned commit.

Distinguish explicitly between:

- `OBSERVED IN SOURCE`;
- `DOCUMENTED`;
- `INFERRED`;
- `UNRESOLVED`.

The MFC 5.0 paper is a capability map, not proof that a combined configuration
is executable.

Never invent command output, successful tests, source locations, or numerical
results.

## Repository boundaries

- Solver-source changes belong in this repository.
- Governing documents, audit reports, reference physics, validation studies, and
  provenance belong in `../MFCNucleationCavitation`.
- Independent C++ oracle work belongs in `../SixEquationCavitation`.
- Do not copy implementation code among the repositories merely to avoid
  understanding an interface.
- Do not modify either sibling repository unless the user's task explicitly
  includes it.

## Git safety

- Do not commit, push, merge, rebase, tag, or open a pull request unless
  explicitly requested.
- Do not modify `origin` or `upstream`.
- Do not force-push.
- Do not reset or rewrite `master`.
- Do not synchronize with a newer upstream commit during an active pinned audit
  unless the user explicitly repins the project baseline.
- Preserve the ability to reproduce commit
  `ca26c4cd48c3d98927a2879a3fc16e53e546a555`.

This `AGENTS.md` is a fork-local research guardrail. Do not include it in an
upstream MFC pull request unless the user and MFC maintainers explicitly decide
otherwise.

## Completion report

After every MFC task, report:

- current branch and commit;
- files inspected;
- files modified;
- commands run;
- tests run and their actual results;
- observed evidence;
- unresolved questions;
- `git status --short`;
- the smallest justified next action.
