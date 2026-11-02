"""Controls for runtime-metadata normalization; no CFD references required."""
import pathlib
import tempfile
import unittest
import xml.etree.ElementTree as ET

import compare_slow_assertions as comparison


class ComparisonTest(unittest.TestCase):
    def compare(self, before, after, expanded=("1 == 1", "1 == 1"), passed=(True, True)):
        with tempfile.TemporaryDirectory() as directory:
            results = []
            for i, contexts in enumerate((before, after)):
                root = ET.Element("doctest")
                case = ET.SubElement(root, "TestCase", name="CLI parity")
                expression = ET.SubElement(case, "Expression", filename="tests/example.cpp", line="1",
                                           type="CHECK", success=str(passed[i]).lower())
                ET.SubElement(expression, "Original").text = "actual == expected"
                ET.SubElement(expression, "Expanded").text = expanded[i]
                for context in contexts:
                    ET.SubElement(expression, "Info").text = context
                ET.SubElement(case, "OverallResultsAsserts", successes=str(int(passed[i])),
                              failures=str(int(not passed[i])))
                path = pathlib.Path(directory) / str(i)
                ET.ElementTree(root).write(path)
                results.append(comparison.read_assertions(path))
            return comparison.compare(*results)

    def test_runtime_metadata_only(self):
        before = """iter continuity x-velocity y-velocity z-velocity k omega sec/iter
1 1e-2 2e-3 3e-4 4e-5 5e-6 6e-7 1.2e-1
Wrote /tmp/kutta-mesh-123-1/fluent (1 iterations, 0.12 s solve)."""
        after = before.replace("1.2e-1", "9.8e-1").replace("123-1", "456-2").replace("0.12", "0.98")
        result = self.compare(["0x101010101", before], ["0x202020202", after])
        self.assertTrue(result["no_new_failures"])
        self.assertEqual(result["changed_values"], [])

    def test_residual_column_remains_part_of_identity(self):
        before = "iter continuity x-velocity sec/iter\n1 1e-2 2e-3 0.1"
        result = self.compare([before], [before.replace("2e-3", "2.00001e-3")])
        self.assertFalse(result["no_new_failures"])
        self.assertEqual(len(result["missing_assertions"]), 1)

    def test_iteration_remains_part_of_identity(self):
        before = "iter continuity x-velocity sec/iter\n1 1e-2 2e-3 0.1"
        self.assertFalse(self.compare([before], [before.replace("\n1 ", "\n2 ")])["no_new_failures"])

    def test_total_iteration_count_is_preserved(self):
        self.assertFalse(self.compare(["Wrote out (3 iterations, 0.1 s solve)"],
                                      ["Wrote out (4 iterations, 0.2 s solve)"])["no_new_failures"])

    def test_model_context_is_preserved(self):
        self.assertFalse(self.compare(["rotating-wheel"], ["pressure-inlet"])["no_new_failures"])

    def test_pointer_aliases_are_preserved(self):
        self.assertFalse(self.compare(["0x101010101", "0x202020202"],
                                      ["0x303030303", "0x303030303"])["no_new_failures"])

    def test_expanded_values_are_compared_exactly(self):
        result = self.compare([], [], expanded=("0.125 == 0.125", "0.125001 == 0.125"))
        self.assertEqual(len(result["changed_values"]), 1)

    def test_changed_outcome_is_a_regression(self):
        self.assertEqual(len(self.compare([], [], passed=(True, False))["regressions"]), 1)

    def test_other_tables_are_not_normalized(self):
        self.assertFalse(self.compare(["iter pressure temperature\n1 2 3"],
                                      ["iter pressure temperature\n1 2 4"])["no_new_failures"])


if __name__ == "__main__":
    unittest.main()
