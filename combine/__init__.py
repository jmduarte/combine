import importlib
import os
from ctypes import cdll

spec = importlib.util.find_spec("_combine")
library_file = os.path.join(os.path.dirname(spec.origin), "libCombine.so")

cdll.LoadLibrary(library_file)

from _combine import _combine as cpp_combine


def combine(
    datacard,
    name="Test",
    method="AsymptoticLimits",
    hint_method="",
    dataset="data_obs",
    toys_file="",
    significance=False,
    verbose=0,
    mass=120,
    systematics=1,
    toys=0,
    expect_signal=0,
    seed=123456,
    lower_limit=False,
    bypass_frequentist_fit=False,
    perf_counters=False,
):
    out = cpp_combine(
        datacard,
        name,
        dataset,
        method,
        hint_method,
        toys_file,
        verbose,
        significance,
        mass,
        systematics,
        toys,
        expect_signal,
        seed,
        lower_limit,
        bypass_frequentist_fit,
        perf_counters,
    )

    return out

def significance(datacard, pvalue=False):
    import scipy.stats
    significance = combine(datacard, method="Significance", mass=125, toys=-1, expect_signal=1, significance=True).getLimit()[-1]

    if pvalue:
        return scipy.stats.norm.sf(significance)

    return significance
