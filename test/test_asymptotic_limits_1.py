import unittest
import pandas as pd
import numpy as np
import combine


class TestAsymptoticLimits1(unittest.TestCase):
    def test_asymptotic_limits_1_01(self):
        """combine -M AsymptiticLimits datacards/simple-counting-experiment.txt
        """
        ref = pd.read_csv("test/outputs/output01.csv", comment="#")

        res = combine.asymptotic_limits("test/datacards/simple-counting-experiment.txt")
        res = res.fillna(0.0)

        np.testing.assert_array_almost_equal(ref["limit"].values, res["limit"].values)
        np.testing.assert_array_almost_equal(ref["limitErr"].values, res["uncertainty"].values)
        np.testing.assert_array_almost_equal(ref["quantileExpected"].values, res["quantile"].values)

    def test_asymptotic_limits_1_02(self):
        """combine -n LimitTest -m 120 -M AsymptoticLimits hgg_8TeV_MVA_cat0145.txt --run both
        """
        ref = pd.read_csv("test/outputs/output02.csv", comment="#")

        res = combine.asymptotic_limits("test/datacards/hgg_8TeV_MVA_cat0145.txt", mass=120)
        res = res.fillna(0.0)

        np.testing.assert_array_almost_equal(ref["limit"].values, res["limit"].values)
        np.testing.assert_array_almost_equal(ref["limitErr"].values, res["uncertainty"].values)
        np.testing.assert_array_almost_equal(ref["quantileExpected"].values, res["quantile"].values)


if __name__ == "__main__":

    unittest.main(verbosity=2)
