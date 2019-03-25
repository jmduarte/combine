import ROOT
import numpy as np

from combine.algos import binned_shape_significance

from sys import exit
from optparse import OptionParser

import tempfile
import os

a = np.random.normal(size=100) + 5
b = np.random.uniform(0, 10, size=1000)

processes = dict(signal=a, background=b)

integrals = dict(signal=10.0, background=100.0)

print(binned_shape_significance(processes, "signal", 0, 10, 10, integrals=integrals))
