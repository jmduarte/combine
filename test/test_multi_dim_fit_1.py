import unittest
import pandas as pd
import numpy as np
import combine


class TestMultiDimFit1(unittest.TestCase):
    def test_multi_dim_fit_1_01(self):
        """-M MultiDimFit -m 125 hgg_8TeV_MVA_cat0145.txt --algo=grid --points 30 --setParameterRanges r=0.0,3.0
        """
        ref = pd.read_csv("test/outputs/output05.csv", comment="#")

        r, delta_nll = combine.multi_dim_fit(
            "test/datacards/hgg_8TeV_MVA_cat0145.txt",
            mass=125,
            algo="grid",
            points=30,
            parameter_ranges="r=0.0,3.0",
        )

        np.testing.assert_array_almost_equal(ref["r"].values, r, decimal=4)
        np.testing.assert_array_almost_equal(ref["deltaNLL"].values, delta_nll, decimal=4)


if __name__ == "__main__":

    unittest.main(verbosity=2)
