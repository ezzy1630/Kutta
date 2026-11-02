#pragma once

// How the app shows file paths. A path on screen ends up in screenshots,
// screen shares and bug reports, so the user's home folder is always written
// "~". Narrow fields show the compact form and keep the full path for a tooltip.

#include <QString>

namespace kapp {

// "/Users/user/cases/car.kutta.json" -> "~/cases/car.kutta.json", with native
// separators. Paths outside home are only made native.
QString homeRelativePath(const QString& path, const QString& home = {});

// The home-relative path with leading folders folded into "…" so that at most
// keepDirs folders are left before the file name:
// "~/…/fsae/cases/car.kutta.json". The file name is never shortened.
QString compactPath(const QString& path, int keepDirs = 2, const QString& home = {});

// Free text (console lines, CLI messages) with every path under home written
// from "~". A home of "/" or "" leaves the text unchanged.
QString withHomeAsTilde(const QString& text, const QString& home = {});

} // namespace kapp
