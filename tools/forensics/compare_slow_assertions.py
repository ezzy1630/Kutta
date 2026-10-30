#!/usr/bin/env python3
"""Compare doctest XML with --success, preserving each assertion's identity.

Usage: compare_slow_assertions.py baseline.xml final.xml > comparison.json
A repeated assertion is identified by its occurrence within the same case,
source location, original expression and INFO context. Changed expanded values
are recorded even when the assertion outcome is unchanged; timings are not compared.
Pointer addresses in INFO are numbered by first appearance within each case;
their process-specific addresses are not stable assertion identities.
Missing assertions and newly failed assertions return a nonzero exit status.
"""
import argparse
import collections
import json
import re
import xml.etree.ElementTree as ET


def read_assertions(path):
    result = {}
    for case in ET.parse(path).getroot().iter("TestCase"):
        if case.get("skipped") == "true":
            continue
        if case.get("name") in result:
            raise ValueError(f"{path}: duplicate case name {case.get('name')}")
        assertions = {}
        occurrences = collections.Counter()
        pointers = {}
        def context_text(info):
            value = " ".join(" ".join(info.itertext()).split())
            def address(match):
                return pointers.setdefault(match.group(), f"<pointer-{len(pointers)+1}>")
            if case.get("name", "").startswith("V7 SA TMR") and re.fullmatch(r"case 0x[0-9a-fA-F]{9,16} level [0-9]+", value):
                return re.sub(r"0x[0-9a-fA-F]{9,16}", address, value)
            return value

        for expression in case.iter("Expression"):
            filename = expression.get("filename", "").replace("\\", "/")
            if "/tests/" in filename:
                filename = "tests/" + filename.split("/tests/", 1)[1]
            context = tuple(context_text(info) for info in expression.findall("Info"))
            key = (filename, expression.get("line"), expression.get("type"),
                   (expression.findtext("Original") or "").strip(), context)
            identity = (*key, occurrences[key])
            occurrences[key] += 1
            assertions[identity] = {
                "passed": expression.get("success") == "true",
                "expanded": (expression.findtext("Expanded") or expression.findtext("Exception") or "").strip(),
            }
        # Exceptions inside a CHECK_NOTHROW/THROWS expression are assertion
        # outcomes. Doctest marks unhandled test-case exceptions with `crash`.
        if any("crash" in error.attrib for error in case.iter("Exception")):
            raise ValueError(f"{path}: interrupted case {case.get('name')}")
        summary = case.find(".//OverallResultsAsserts")
        if summary is None:
            raise ValueError(f"{path}: incomplete case {case.get('name')}")
        expected = int(summary.get("successes", 0)) + int(summary.get("failures", 0))
        if len(assertions) != expected:
            raise ValueError(f"{path}: {case.get('name')} lacks assertions; run with --success")
        result[case.get("name")] = assertions
    if not result:
        raise ValueError(f"{path}: no executed test cases")
    return result


def compare(before, after):
    rows, regressions, recoveries, missing, added, changed = [], [], [], [], [], []
    for name, assertions in before.items():
        final = after.get(name, {})
        rows.append({"case": name, "baseline_assertions": len(assertions),
                     "final_assertions": len(final),
                     "baseline_failed": sum(not a["passed"] for a in assertions.values()),
                     "final_failed": sum(not a["passed"] for a in final.values())})
        for key, old in assertions.items():
            if key not in final:
                missing.append({"case": name, "assertion": key})
                continue
            new = final[key]
            if old["expanded"] != new["expanded"]:
                changed.append({"case": name, "assertion": key,
                                "before": old["expanded"], "after": new["expanded"],
                                "baseline_passed": old["passed"], "final_passed": new["passed"]})
            if old["passed"] != new["passed"]:
                change = {"case": name, "assertion": key,
                          "before": old["expanded"], "after": new["expanded"]}
                (regressions if old["passed"] else recoveries).append(change)
    for name, assertions in after.items():
        for key, value in assertions.items():
            if key not in before.get(name, {}):
                added.append({"case": name, "assertion": key, **value})
                if not value["passed"]:
                    regressions.append({"case": name, "assertion": key,
                                        "before": None, "after": value["expanded"]})
    missing_cases = sorted(before.keys() - after.keys())
    return {"cases": rows, "regressions": regressions, "recoveries": recoveries,
            "changed_values": changed,
            "missing_assertions": missing, "new_assertions": added,
            "missing_cases": missing_cases,
            "new_cases": sorted(after.keys() - before.keys()),
            "no_new_failures": not regressions and not missing and not missing_cases}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline")
    parser.add_argument("final")
    args = parser.parse_args()
    try:
        result = compare(read_assertions(args.baseline), read_assertions(args.final))
    except (OSError, ET.ParseError, ValueError) as error:
        parser.error(str(error))
    print(json.dumps(result, indent=2))
    return 0 if result["no_new_failures"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
