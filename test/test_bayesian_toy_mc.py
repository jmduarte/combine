import unittest
import pandas as pd
import numpy as np
import combine


class TestBayesianToyMC(unittest.TestCase):
    def test_bayesian_toy_mc_01(self):
        """combine -M BayesianToyMC datacards/simple-counting-experiment.txt
        """
        res = combine.bayesian_toy_mc("test/datacards/simple-counting-experiment.txt")

        self.assertAlmostEqual(res, 0.6734255553838588, places=6)


if __name__ == "__main__":

    unittest.main(verbosity=2)
