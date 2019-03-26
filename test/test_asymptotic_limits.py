import unittest
import pandas as pd
import numpy as np
import combine


class TestAsymptoticLimits(unittest.TestCase):
    def test_asymptotic_limits_01(self):
        """combine -M AsymptiticLimits datacards/simple-counting-experiment.txt
        """
        ref = pd.read_csv("test/outputs/output01.csv", comment="#")

        res = combine.asymptotic_limits("test/datacards/simple-counting-experiment.txt")
        res = res.fillna(0.0)

        np.testing.assert_array_almost_equal(ref["limit"].values, res["limit"].values)
        np.testing.assert_array_almost_equal(ref["limitErr"].values, res["uncertainty"].values)
        np.testing.assert_array_almost_equal(ref["quantileExpected"].values, res["quantile"].values)


if __name__ == "__main__":

    unittest.main(verbosity=2)
