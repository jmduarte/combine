import unittest
import pandas as pd
import numpy as np
import combine


class TestSignificance(unittest.TestCase):
    def test_significance_01(self):
        """combine -M Significance datacards/simple-counting-experiment.txt
        """
        res = combine.significance("test/datacards/simple-counting-experiment.txt")

        self.assertEqual(res, 0.0)

    def test_significance_02(self):
        """combine -M Significance datacards/simple-counting-experiment.txt --pval
        """
        res = combine.significance("test/datacards/simple-counting-experiment.txt", pvalue=True)

        self.assertEqual(res, 0.5)

    def test_significance_03(self):
        """combine -M Significance datacards/simple-counting-experiment.txt -t -1 --significance --expectSignal=1
        """
        res = combine.significance(
            "test/datacards/simple-counting-experiment.txt", toys=-1, expect_signal=1, significance=True
        )

        self.assertAlmostEqual(res, 2.36714142758297)

    def test_significance_04(self):
        """combine -M Significance datacards/simple-counting-experiment.txt -t -1 --significance --expectSignal=1 --pval
        """
        res = combine.significance(
            "test/datacards/simple-counting-experiment.txt", toys=-1, expect_signal=1, significance=True, pvalue=True
        )

        self.assertAlmostEqual(res, 0.008963040597835892)


if __name__ == "__main__":

    unittest.main(verbosity=2)
