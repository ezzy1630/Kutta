#pragma once

#include <QString>

namespace kapp {

// Match score of a typed query against a label: prefix > substring >
// subsequence (word starts and runs score higher). < 0 means no match.
int fuzzyScore(const QString& query, const QString& text);

} // namespace kapp
