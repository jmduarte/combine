import unittest
import pandas as pd
import numpy as np
import combine


class TestFitDiagnostics1(unittest.TestCase):
    def test_fit_diagnostics_1_01(self):
        """-M FitDiagnostics -m 125 hgg_8TeV_MVA_cat0145.txt --robustFit=1
        """
        ref = pd.read_csv("test/outputs/output04.csv", comment="#")

        res = combine.fit_diagnostics("test/datacards/hgg_8TeV_MVA_cat0145.txt", mass=125, robust_fit=True)
        res = res.fillna(0.0)

        np.testing.assert_array_almost_equal(ref["limit"].values, res["limit"].values)
        np.testing.assert_array_almost_equal(ref["limitErr"].values, res["uncertainty"].values)
        np.testing.assert_array_almost_equal(ref["quantileExpected"].values, res["quantile"].values)


if __name__ == "__main__":

    unittest.main(verbosity=2)
