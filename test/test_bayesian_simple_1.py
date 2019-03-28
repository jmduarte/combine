import unittest
import pandas as pd
import numpy as np
import combine


class TestBayesianSimple(unittest.TestCase):
    def test_bayesian_simple_01(self):
        """combine -M BayesianSimple datacards/simple-counting-experiment.txt
        """
        res = combine.bayesian_simple("test/datacards/simple-counting-experiment.txt")

        self.assertAlmostEqual(res, 0.6722917026656567, places=6)


if __name__ == "__main__":

    unittest.main(verbosity=2)
