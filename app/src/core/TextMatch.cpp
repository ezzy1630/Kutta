#include "core/TextMatch.h"

#include <algorithm>

namespace kapp {

int fuzzyScore(const QString& query, const QString& text)
{
    if (query.isEmpty()) return 0;
    const QString q = query.toLower(), t = text.toLower();
    // Whole-word prefix hits rank first, then contiguous substrings, then subsequences.
    if (t.startsWith(q)) return 1000 - int(t.size());
    if (int i = int(t.indexOf(q)); i >= 0) return 800 - i - int(t.size()) / 4;
    int score = 0, ti = 0, run = 0;
    for (const QChar c : q) {
        if (c == ' ') continue;
        const int found = int(t.indexOf(c, ti));
        if (found < 0) return -1;
        run = found == ti ? run + 1 : 0;
        const bool boundary = found == 0 || !t[found - 1].isLetterOrNumber();
        score += 10 + run * 5 + (boundary ? 15 : 0) - std::min(9, found - ti);
        ti = found + 1;
    }
    return score;
}

} // namespace kapp
