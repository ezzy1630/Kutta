# Desktop app dependencies and licences

Licences do not drive choices for now (one engineer's use, not a product). This
list exists so a later commercial decision starts from facts. One line each;
"bundled" means it ships inside the packaged app.

- Qt 6.12 Core, Gui, Widgets, Svg (bundled): LGPL-3.0, GPL or commercial; dynamically linked.
- Qt GuiPrivate / QRhi (bundled with QtGui): same licence; semi-public API, so the app ships with the exact Qt it was built against.
- Qt ShaderTools (build only, `qsb` shader compilation): GPL-3.0 or commercial; not shipped.
- Qt Test (tests only): LGPL-3.0, GPL or commercial; not shipped.
- Qt Advanced Docking System 5.1.1 (bundled, shared library): LGPL-2.1; text in the package's `licenses/`.
- Lucide icons (bundled, compiled into the app): ISC; `app/resources/icons/LICENSE-lucide.txt`.
- Colormap samples: viridis, cividis, inferno from matplotlib (matplotlib licence, BSD-style); cool–warm after Moreland (public description); RdBu from ColorBrewer (Apache-2.0).
- Ansys mixing elbow example (recorded fixture data in `app/fixtures/`, test data in `app/tests/data/`): MIT; see `validation/LICENSE.ansys-example-data`.
- kutta CLI (bundled as a separate executable): this repository's own code.
