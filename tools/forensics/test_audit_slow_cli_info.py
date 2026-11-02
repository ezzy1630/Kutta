import unittest

from audit_slow_cli_info import CASES, normalize, normalize_info
from compare_slow_assertions import compare


class CliInfoAuditTests(unittest.TestCase):
    def log(self, seconds, residual='1.000000e-03'):
        rows = ' '.join(f'{i} {residual} 2.000000e-03 {seconds}' for i in (1, 2, 3))
        return ('iter continuity x-velocity sec/iter ' + rows
                + ' Reached the iteration limit (3 iterations). Wrote '
                + '/var/folders/xx/task/T/kutta-mesh-1234-1/fluent '
                + f'(3 iterations, {seconds} s solve).')

    def test_clock_normalization_preserves_physical_residuals(self):
        a = normalize_info(self.log('1.000000e-01'), {})
        b = normalize_info(self.log('9.000000e-01'), {})
        self.assertEqual(a, b)
        self.assertNotEqual(a, normalize_info(self.log('9.000000e-01', '1.100000e-03'), {}))
        with self.assertRaises(ValueError):
            normalize_info(self.log('1.000000e-01').replace('Reached the iteration limit', 'Unexpected exit'), {})

    def test_pointer_aliasing_is_preserved(self):
        a, b = {}, {}
        first = [normalize_info(p, a) for p in ('0x100001111', '0x100002222', '0x100001111')]
        second = [normalize_info(p, b) for p in ('0x200001111', '0x200002222', '0x200002222')]
        self.assertNotEqual(first, second)

    def test_expanded_values_and_failures_are_never_normalized(self):
        case = sorted(CASES)[0]
        key = ('tests/example.cpp', '10', 'CHECK', 'physical_error < tolerance', ('0x100001111',), 0)
        baseline = {case: {key: {'passed': True, 'expanded': '0.09 < 0.10'}}}
        branch_key = (*key[:4], ('0x200001111',), 0)
        branch = {case: {branch_key: {'passed': False, 'expanded': '0.11 < 0.10'}}}
        result = compare(normalize(baseline), normalize(branch))
        self.assertEqual(len(result['regressions']), 1)
        self.assertEqual(len(result['changed_values']), 1)
        self.assertEqual(len(result['missing_assertions']), 0)


if __name__ == '__main__':
    unittest.main()
