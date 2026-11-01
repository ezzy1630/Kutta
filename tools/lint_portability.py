#!/usr/bin/env python3
"""Fail on source patterns known to break one of Kutta's supported toolchains.

Floating-point std::from_chars is missing from Xcode 16's libc++ (macOS CI) but
compiles on newer local toolchains, so it slips through local builds. Use
kutta::parse_double instead. Integer from_chars is fine everywhere.

Absolute home-directory paths only resolve on the machine that wrote them, so
fixture tests elsewhere silently skip. Build paths from KUTTA_SOURCE_DIR. In any
tracked text (docs, logs, results) they also publish the author's local layout
from this public repository; write <checkout> instead.

std::jthread is likewise missing from that libc++; use std::thread and join.

Tests reach the Fluent corpus only through tests/fluent_examples.hpp, which
honours KUTTA_FLUENT_EXAMPLES; a hard-coded corpus path silently skips in any
checkout (worktree, CI) that keeps the corpus elsewhere.

HDF5 takes char* file names and expects UTF-8; path::c_str() is wchar_t on
Windows and path::string() uses the ANSI code page there. Use kutta::utf8_path.
"""

import pathlib
import re
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
SOURCE_DIRS = ("src", "tests", "bench", "tools")
SUFFIXES = {".cpp", ".hpp", ".h", ".cu", ".mm"}
# Third argument of from_chars(first, last, value[, base_or_format]).
CALL = re.compile(r"from_chars\s*\([^;]*?,[^;]*?,\s*(\w+)\s*[,)]")
JTHREAD = re.compile(r"\bstd::jthread\b")
CORPUS = re.compile(r"validation/fluent_examples\b")
H5_FILE = re.compile(r"\bH5F(?:open|create|is_accessible)\s*\(\s*(?!kutta::utf8_path\()")
# Generic placeholder accounts (upstream docs, CI runners) reveal nothing.
ANY_HOME = re.compile(r"(?:/Users|/home)/(?!(?:user|username|runner)/)[A-Za-z0-9._-]+/"
                      r"|[A-Za-z]:\\{1,2}Users\\{1,2}(?!(?:user|username|runneradmin)\\)")
HOME_PATH = re.compile(r"[\"']((?:/Users|/home)/[^\"']*|[A-Za-z]:\\\\Users\\\\[^\"']*)")


def float_declared(text: str, name: str, before: int) -> bool:
    # The nearest declaration of `name` before the call decides its type.
    decls = list(re.finditer(r"\b([A-Za-z_][\w:<>]*)\s+" + re.escape(name) + r"\s*[=;{(,]", text[:before]))
    decls = [d for d in decls if d.group(1) not in {"return", "auto&", "const"}]
    return bool(decls) and decls[-1].group(1) in {"double", "float"}


def tracked_files() -> list[str]:
    tracked = subprocess.run(["git", "-C", str(ROOT), "ls-files", "-z"], capture_output=True)
    if tracked.returncode == 0:
        return tracked.stdout.decode().split("\0")
    # CI containers without git check out a plain tree, which holds only tracked files.
    return [str(p.relative_to(ROOT)) for p in ROOT.rglob("*")
            if not any(part == ".git" or part.startswith("build") for part in p.relative_to(ROOT).parts)]


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
            for m in H5_FILE.finditer(text) if path.suffix != ".py" else ():
                line = text.count("\n", 0, m.start()) + 1
                bad.append(f"{path.relative_to(ROOT)}:{line}: HDF5 file name must be kutta::utf8_path(...)")
            for m in JTHREAD.finditer(text) if path.suffix != ".py" else ():
                line = text.count("\n", 0, m.start()) + 1
                bad.append(f"{path.relative_to(ROOT)}:{line}: std::jthread is missing from macOS CI libc++; use std::thread")
            if d == "tests" and path.suffix != ".py" and path.name != "fluent_examples.hpp":
                for m in CORPUS.finditer(text):
                    line = text.count("\n", 0, m.start()) + 1
                    bad.append(f"{path.relative_to(ROOT)}:{line}: Fluent corpus path; use fluent_examples() "
                               f"from tests/fluent_examples.hpp")
            if path.suffix == ".py" or path.name == "parse.cpp":
                continue
            for m in CALL.finditer(text):
                if float_declared(text, m.group(1), m.start()):
                    line = text.count("\n", 0, m.start()) + 1
                    bad.append(f"{path.relative_to(ROOT)}:{line}: floating-point std::from_chars "
                               f"into '{m.group(1)}'; use kutta::parse_double")
    for name in tracked_files():
        path = ROOT / name
        if not name or not path.is_file():
            continue
        data = path.read_bytes()
        if b"\0" in data[:8192]:
            continue  # binary
        text = data.decode(errors="replace")
        for m in ANY_HOME.finditer(text):
            line = text.count("\n", 0, m.start()) + 1
            bad.append(f"{name}:{line}: local home path '{m.group(0)}'; write <checkout> instead")
    for b in bad:
        print(b, file=sys.stderr)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
