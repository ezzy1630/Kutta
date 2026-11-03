"""Budget controls: physical omega destruction and reference dump identity."""
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

import sst_budget


class BudgetTests(unittest.TestCase):
    def test_newton_compensation_and_signed_cross_diffusion(self):
        names = 'cell k omega f1 cdkw pk dk ps ds ck dk_transport cw dw rk rw fixed'.split()
        cells = np.zeros(2, dtype=[(name, float) for name in names])
        cells['cell'] = [0, 1]
        cells['k'], cells['omega'], cells['f1'] = 1, 2, .4
        cells['pk'], cells['dk'] = .18, .18
        cross = np.array([.6, -.7])
        production = np.array([1., 2.])
        destruction = (.075 * .4 + .0828 * .6) * 4
        cells['cdkw'] = cross / .6
        cells['ps'] = production + np.maximum(cross, 0) + destruction
        cells['ds'] = (2 * destruction + np.maximum(-cross, 0)) / 2
        cells['rw'] = production + cross - destruction
        st = SimpleNamespace(rho=np.ones(2), ctx=SimpleNamespace(vol=np.ones(2)))
        result = sst_budget.budget(cells, st, 'default')['omega/interior']
        self.assertAlmostEqual(result['Y'], 2 * destruction)
        self.assertAlmostEqual(result['P/Y'], production.sum() / (2 * destruction))
        self.assertAlmostEqual(result['X/Y'], cross.sum() / (2 * destruction))
        self.assertAlmostEqual(result['split_roundoff/Y'], 0)
        self.assertAlmostEqual(result['native_absR/Y'], result['absR/Y'])

    def test_cell_identity_reorders_and_rejects_invalid_dumps(self):
        with tempfile.TemporaryDirectory() as root:
            directory = Path(root)
            identity = 'case\ndata\ncommit\nmetal\nsaved\ndefault\n0\n'
            (directory / 'identity.txt').write_text(identity)
            (directory / 'cells.csv').write_text('cell,k\n1,2\n0,1\n')
            cells, _ = sst_budget.native_cells(directory, ['case', 'data'], 2)
            np.testing.assert_equal(cells['k'], [1, 2])
            with self.assertRaisesRegex(ValueError, 'identity'):
                sst_budget.native_cells(directory, ['other-case', 'data'], 2)
            for content, reason in [('cell,k\n0,1\n0,2\n', 'exactly once'),
                                    ('cell,k\n0,1\n1,nan\n', 'nonfinite')]:
                (directory / 'cells.csv').write_text(content)
                with self.assertRaisesRegex(ValueError, reason):
                    sst_budget.native_cells(directory, ['case', 'data'], 2)


if __name__ == '__main__':
    unittest.main()
