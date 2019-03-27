import unittest
import pandas as pd
import numpy as np
import combine


class TestSignificance(unittest.TestCase):
    def test_significance_07(self):
        """-M Significance --significance -m 125 hgg_8TeV_MVA_cat0145.txt -t -1 --expectSignal=1 --toysFreq
        """
        res = combine.significance(
            "test/datacards/hgg_8TeV_MVA_cat0145.txt",
            toys=-1,
            expect_signal=1,
            significance=True,
            mass=125,
            toys_frequentist=True,
        )

        self.assertAlmostEqual(res, 2.7649763354677384, places=6)


if __name__ == "__main__":

    unittest.main(verbosity=2)
