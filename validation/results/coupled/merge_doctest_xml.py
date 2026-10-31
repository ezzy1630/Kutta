"""Merge doctest XML reports of disjoint test-case subsets into one report.

Usage: python3 merge_doctest_xml.py out.xml.gz part1.xml[.gz] part2.xml[.gz] ...

The slow tier was run in chunks (-tc / -tce filters) so that no single run
outlives a background-job limit; compare_slow_assertions.py reads one file.
A case reported by two parts is an error.
"""
import gzip
import sys
import xml.etree.ElementTree as ET


def load(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rb") as f:
        return ET.parse(f).getroot()


out, parts = sys.argv[1], sys.argv[2:]
root = load(parts[0])
group = root.find("TestSuite")
seen = {case.get("name") for case in root.iter("TestCase") if case.get("skipped") != "true"}
for path in parts[1:]:
    for suite in load(path).iter("TestSuite"):
        for case in suite.findall("TestCase"):
            if case.get("skipped") == "true":
                continue
            name = case.get("name")
            if name in seen:
                raise SystemExit(f"{path}: case {name!r} already present")
            seen.add(name)
            group.append(case)  # skipped placeholders are ignored by the comparison
data = ET.tostring(root, encoding="utf-8", xml_declaration=True)
with (gzip.open if out.endswith(".gz") else open)(out, "wb") as f:
    f.write(data)
print(f"{len(seen)} cases")
