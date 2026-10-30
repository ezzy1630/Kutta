---
version: 1
slug: "src-main-cpp"
primary_target: "src/main.cpp"
related_targets: []
---

# Kutta desktop workbench: surface brief

Scope: the whole native workbench (outline, inspector, graphics, charts, console, palette, history). Mode: Operate.
Audience: experienced Fluent user (FSAE aero engineer) on Windows/RTX; task: run their own Fluent case and decide whether Kutta's forces and fields agree with Fluent.
Constraints: Fluent-recognisable structure and vocabulary (owner-pinned); no idle repaint; approximate/modified/synthetic labelled everywhere; light and dark; crisp at 100/125/150%.

## Direction contract

THESIS: A race engineer's data-analysis workbench wearing Fluent's outline. Refuses the stock grey Qt/SaaS shell and Fluent's ribbon-and-dialog clutter: one window, docked panes, every value inline.

OWN-WORLD: Telemetry software (MoTeC i2 / AiM) as the hand: graphite panels in dark, paper-white in light, hairline separators, no cards or shadows, tabular figures in the system UI face. Chrome is neutral; colour lives only in data plates (viewports, charts, colour bars) and state. One instrument accent, light-sheet green, marks the cursor, focus, selection and active probe. Sodium amber is reserved for APPROXIMATE/MODIFIED, red for blocked/error.

STORY: Open the case, see its verdict and zones in Fluent's tree, change a boundary value inline, run, watch residuals fall, then put Kutta beside Fluent and read the delta.

FIRST VIEWPORT: Toolbar with case, verdict badge, device, Run/Stop and iteration field. Left: outline tree. Centre: graphics viewport(s) owning the space. Right: inspector. Bottom: residuals chart and console, tabbed.

FORM: seed: bda2b260d2c26895da2271d9403e6163, assigned grounded candidate 3 (motorsport telemetry analysis). Raises: changed values hold their marker until reset (from gate-board); dense ruled tables (split-flap board).
Signature move: the linked cursor. Hovering a chart sample, XY-plot point or 3D pick shows one green crosshair and readout in every view showing the same quantity.

FINISH: screenshots of every major screen in both themes in docs/app-screens/, each carrying its provenance (fixture vs real CLI, synthetic labels).
