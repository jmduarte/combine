import unittest
import pandas as pd
import numpy as np
import combine


class TestSignificance(unittest.TestCase):
    def test_significance_05(self):
        """-M significance --significance -m 125 hgg_8TeV_MVA_cat0145.txt
        """
        res = combine.significance("test/datacards/hgg_8TeV_MVA_cat0145.txt", significance=True, mass=125)

        self.assertAlmostEqual(res, 2.246109366413282)


if __name__ == "__main__":

    unittest.main(verbosity=2)
