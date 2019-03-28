import unittest
import pandas as pd
import numpy as np
import combine


class TestHybridNew(unittest.TestCase):
    def test_hybrid_new_01(self):
        """combine -M HybridNew datacards/simple-counting-experiment.txt
        """
        limit, error = combine.hybrid_new("test/datacards/simple-counting-experiment.txt")

        self.assertAlmostEqual(limit, 0.640585647004483, places=6)
        self.assertAlmostEqual(error, 0.013576435560588898, places=6)


if __name__ == "__main__":

    unittest.main(verbosity=2)
