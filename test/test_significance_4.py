import unittest
import pandas as pd
import numpy as np
import combine


class TestSignificance(unittest.TestCase):
    def test_significance_06(self):
        """-M Significance --significance -m 125 hgg_8TeV_MVA_cat0145.txt -t -1 --expectSignal=1
        """
        res = combine.significance(
            "test/datacards/hgg_8TeV_MVA_cat0145.txt", toys=-1, expect_signal=1, significance=True, mass=125
        )

        self.assertAlmostEqual(res, 2.729286751814355, places=6)


if __name__ == "__main__":

    unittest.main(verbosity=2)
