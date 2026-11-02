#pragma once

// Checks scripted sessions (--scenario) run on what is on screen before a
// screenshot, so the smoke test and every capture catch what a reviewer would:
// text that overlaps, is clipped or is cut short, and the user's home folder.

#include <QStringList>

class QWidget;

namespace kapp {

// Geometry problems among root's visible text and controls, as laid out now:
// a control narrower than its minimum, a wrapped label shorter than its text,
// a control outside the widget that holds it, or two controls overlapping.
QStringList layoutProblems(QWidget* root);

// Visible text under root (labels, buttons, fields, the console, table cells,
// tabs) that contains the home folder. Tooltips are not on screen.
QStringList homePathsOnScreen(QWidget* root);

} // namespace kapp
