#!/usr/bin/env python3
"""Fail on source patterns known to break one of Kutta's supported toolchains.

Floating-point std::from_chars is missing from Xcode 16's libc++ (macOS CI) but
compiles on newer local toolchains, so it slips through local builds. Use
kutta::parse_double instead. Integer from_chars is fine everywhere.

Absolute home-directory paths only resolve on the machine that wrote them, so
fixture tests elsewhere silently skip. Build paths from KUTTA_SOURCE_DIR.
"""

import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
SOURCE_DIRS = ("src", "tests", "bench", "tools")
SUFFIXES = {".cpp", ".hpp", ".h", ".cu", ".mm"}
# Third argument of from_chars(first, last, value[, base_or_format]).
CALL = re.compile(r"from_chars\s*\([^;]*?,[^;]*?,\s*(\w+)\s*[,)]")
HOME_PATH = re.compile(r"[\"']((?:/Users|/home)/[^\"']*|[A-Za-z]:\\\\Users\\\\[^\"']*)")


def float_declared(text: str, name: str, before: int) -> bool:
    # The nearest declaration of `name` before the call decides its type.
    decls = list(re.finditer(r"\b([A-Za-z_][\w:<>]*)\s+" + re.escape(name) + r"\s*[=;{(,]", text[:before]))
    decls = [d for d in decls if d.group(1) not in {"return", "auto&", "const"}]
    return bool(decls) and decls[-1].group(1) in {"double", "float"}


def main() -> int:
    bad = []
    for d in SOURCE_DIRS:
        for path in sorted((ROOT / d).rglob("*")):
            if path.suffix not in SUFFIXES | {".py"}:
                continue
            text = path.read_text(errors="replace")
            for m in HOME_PATH.finditer(text):
                line = text.count("\n", 0, m.start()) + 1
                bad.append(f"{path.relative_to(ROOT)}:{line}: absolute home path '{m.group(1)}'; "
                           f"build it from KUTTA_SOURCE_DIR")
            if path.suffix == ".py" or path.name == "parse.cpp":
                continue
            for m in CALL.finditer(text):
                if float_declared(text, m.group(1), m.start()):
                    line = text.count("\n", 0, m.start()) + 1
                    bad.append(f"{path.relative_to(ROOT)}:{line}: floating-point std::from_chars "
                               f"into '{m.group(1)}'; use kutta::parse_double")
    for b in bad:
        print(b, file=sys.stderr)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
