# Product

<!-- impeccable:product-schema 1 -->

Written by the implementation agent from the owner's written brief
(2026-10-30) and the earlier product brief in commit 4699d0b, without an
interview. Facts marked *(inferred)* are the agent's reading of the brief, not
confirmed answers.

## Platform

desktop

Native desktop app for Windows (primary user) and macOS (development machine).
Follows each OS's conventions for menus, shortcuts, file dialogs and high-DPI.

## Stack

Decided by the owner after comparing Qt, Rust GUIs, Dear ImGui, VTK and raw
platform APIs: Qt 6 (latest stable), C++20, built by the same CMake project as
the solver. Qt Widgets for the dense workbench, QRhiWidget for 3D viewports
(Metal on macOS, D3D on Windows), own charts. Licensing is not a constraint
for now; libraries are chosen on quality, performance and maintainability and
listed in `docs/app-licensing.md`. The solver runs as the `kutta` CLI child
process behind `docs/app-protocol.md`; it is never linked into the GUI.

## Users

An experienced Ansys Fluent user. The first is an FSAE (Formula Student) aero
engineer who runs Fluent external-aero cases on a car, on their own Windows PC
with an NVIDIA RTX GPU. They think in Fluent's vocabulary (velocity-inlet,
pressure-outlet, report definitions, scaled residuals, Cl/Cd, y+) and judge a
solver by whether it reproduces Fluent on the same mesh.

## Product Purpose

Make that user switch: open the same Fluent case, get the same results, much
faster, on their own machine. Success is the engineer trusting Kutta's
downforce and drag on their car because the app showed them, side by side,
that Kutta agrees with their own Fluent run. The owner: "similar to Ansys,
just much better and easier to use, but can still have all that control and
different views and charts and data". UI/UX matters as much as the numerics.

## Positioning

Same-mesh reproduction of the user's own Fluent case on one consumer GPU, with
the comparison against their Fluent data built in. It reads the team's case;
it does not mesh and does not replace Fluent's pre-processing.

## Operating Context

- Inputs: Fluent `.cas.h5`, optionally the `.dat.h5` from the same run.
- Runs last seconds to tens of minutes while the GPU is saturated. The app
  must not compete for the GPU: nothing repaints while idle.
- The user moves between this app, Fluent, Excel and slide decks for design
  reviews; CSV, PNG and VTU exports feed those.
- Windows at 100/125/150% display scaling is the primary target *(inferred)*;
  macOS on Apple Silicon is the development machine.

## Capabilities and Constraints

- Fluent-style workbench: outline tree in Fluent's order and vocabulary
  (Setup, Solution, Results), inspector panels, 1/2/4 graphics viewports,
  residual/report/XY charts, console, command palette, run history and
  comparison. Recognisable to a Fluent user instantly, minus Fluent's clutter,
  modal dialogs and hidden state.
- The CLI is the single source of truth for what a case means. Editable
  settings come from its `settings` inventory and are applied through
  `solve --setup` overrides; anything else is shown read-only with its reason.
- Support verdict per case: supported, approximate (runs only after explicit
  opt-in per approximation), unsupported (blocked). An approximate result must
  never look faithful anywhere: badges, headers, history, exports, screenshots.
- Convergence criteria (Fluent's residual criteria and report convergence
  conditions) are off by default, as in older scripts: a run stops at its
  iteration cap unless the user turns checking on, and the app claims
  convergence only when the CLI reports which conditions were met.
- Aero numbers (coefficients, groups, balance, centre of pressure) come from
  `kutta aero` over the run's wall loads, so changing groups or the vehicle
  never re-solves. Design studies are the app queueing ordinary runs; each
  variant is a run with its own provenance.

## Evidence on Hand

- Real Fluent 22.2 mixing-elbow case and data (Ansys example data, MIT), and
  real Kutta runs of it (`app/fixtures/cases/mixing_elbow`).
- Measured V6 diagnostic numbers in `src/kutta/fv/VERIFICATION.md`.
- No real car case. Car fixtures are synthetic and labelled so in every file;
  no benchmark, speed-up or agreement claim may be invented from them.

## Product Principles

1. Fluent's words, Fluent's units. Never make a Fluent user translate.
2. Honesty over reassurance: approximate, modified, unconverged and synthetic
   states are always visible and never softened.
3. The comparison with the user's own Fluent run is the product's argument.
4. The GPU belongs to the solver: the UI redraws only on input or new data.
5. Full control without clutter: every setting reachable, nothing modal,
   nothing hidden.

## Accessibility & Inclusion

WCAG AA contrast in both themes; full keyboard operation; accessible names on
every control; colour scales and deltas legible with colour-vision deficiency
*(inferred)*.
