import ROOT

from combine.DatacardParser import *
from combine.ModelTools import *
from combine.ShapeTools import *
from combine.PhysicsModel import *

from sys import modules


def dc2workspace(
    DC,
    filename="./",
    stat=False,
    fixpars=False,
    cexpr=False,
    bin=True,
    out=None,
    verbose=0,
    mass=0,
    dataname="data_obs",
    libs=None,
    modelparams=[],
    poisson=0,
    defMorph="shape",
    noBOnly=False,
    noHistFuncWrappers=False,
    noOptimizePdf=False,
    moreOptimizeSimPdf="none",
    doMasks=False,
    useHistPdf="never",
    nuisancesToRescale=[],
    nuisanceFunctions=[],
    nuisanceGroupFunctions=[],
    forceNonSimPdf=False,
    noCheckNorm=False,
    optimizeExistingTemplates=True,
    optimizeBoundNuisances=True,
    optimizeTemplateBins=True,
    physModel="combine.PhysicsModel:defaultModel",
    physOpt=[],
    dumpCard=False,
):

    if dumpCard:
        DC.print_structure()
        exit()

    ## Load tools to build workspace
    MB = None
    if DC.hasShapes:
        MB = ShapeBuilder(
            DC,
            libs,
            forceNonSimPdf,
            optimizeExistingTemplates,
            moreOptimizeSimPdf,
            poisson,
            modelparams,
            doMasks,
            noHistFuncWrappers,
            noCheckNorm,
            useHistPdf,
            optimizeTemplateBins,
            fixpars,
            defMorph,
            mass,
            filename,
            bin,
            out,
            verbose,
            cexpr,
            dataname,
            nuisanceFunctions,
            nuisanceGroupFunctions,
            nuisancesToRescale,
            noOptimizePdf,
            optimizeBoundNuisances,
            noBOnly,
        )
    else:
        MB = CountingModelBuilder(
            DC,
            mass,
            filename,
            bin,
            out,
            verbose,
            cexpr,
            dataname,
            nuisanceFunctions,
            nuisanceGroupFunctions,
            nuisancesToRescale,
            noOptimizePdf,
            optimizeBoundNuisances,
            noBOnly,
        )

    ## Load physics model
    (physModMod, physModName) = physModel.split(":")
    __import__(physModMod)
    mod = modules[physModMod]
    physics = getattr(mod, physModName)
    if mod == None:
        raise RuntimeError("Physics model module %s not found" % physModMod)
    if physics == None or not isinstance(physics, PhysicsModelBase):
        raise RuntimeError(
            "Physics model %s in module %s not found, or not inheriting from PhysicsModelBase"
            % (physModName, physModMod)
        )
    physics.setPhysicsOptions(physOpt)
    ## Attach to the tools, and run
    MB.setPhysics(physics)
    MB.doModel()


def text2workspace(
    fileName,
    noJMax=False,
    allowNoSignal=False,
    allowNoBackground=False,
    evaluateEdits=True,
    nuisancesToExclude=[],
    stat=False,
    verbose=0,
    **kwargs
):

    if fileName.endswith(".gz"):
        import gzip

        f = gzip.open(fileName, "rb")
        fileName = fileName[:-3]
    else:
        f = open(fileName, "r")

    ## Parse text file
    DC = parseCard(
        f, bin, noJMax, allowNoSignal, allowNoBackground, evaluateEdits, nuisancesToExclude, stat, verbose
    )

    f.close()

    return dc2workspace(DC, stat=stat, verbose=verbose, filename=fileName, **kwargs)



def get_module_file(name):

    import sys

    if sys.version_info < (3, 0):
        import imp

        return imp.find_module(name)[1]

    elif sys.version_info < (3, 5):
        import pkgutil

        package = pkgutil.get_loader(name)
        return package.filename

    import importlib

    spec = importlib.util.find_spec(name)
    return spec.origin


import os
from ctypes import cdll

library_file = os.path.join(os.path.dirname(get_module_file("_combine")), "libCombine.so")

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
    libraries=[],
    compiled=False,
    keyword_value=[],
    text2workspace_kwargs={},
    toys_frequentist=False,
    toys_no_systematics=False,
    robust_fit=False,
    algo=0,
    points=50,
    parameter_ranges=""
):
    if toys_no_systematics and toys_frequentist:
        raise ValueError("You can't set toysNoSystematics and toysFrequentist options at the same time")

    import os

    cwd = os.getcwd()

    if not datacard.startswith("/"):
        datacard = os.path.join(cwd, datacard)

    tmp_file = None
    if not datacard.endswith(".root"):
        import tempfile

        tmp_file = os.path.join(tempfile.gettempdir(), next(tempfile._get_candidate_names()) + ".root")
        text2workspace(
            datacard,
            mass=mass,
            dataname=dataset,
            stat=not systematics,
            verbose=verbose - 1,
            noBOnly=method in ["FitDiagnostics", "MultiDimFit"],
            bin=True,
            out=tmp_file,
            libs=libraries,
            cexpr=compiled,
            modelparams=keyword_value,
            **text2workspace_kwargs
        )

        datacard = tmp_file

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
        toys_frequentist,
        toys_no_systematics,
        robust_fit,
        algo,
        points,
        parameter_ranges
    )

    if not tmp_file is None:
        os.remove(tmp_file)

    return out


def asymptotic_limits(datacard, **kwargs):

    import numpy as np
    import pandas as pd

    out = combine(datacard, method="AsymptoticLimits", **kwargs)

    df = pd.DataFrame(
        data=dict(limit=out.getLimit(), uncertainty=out.getLimitError(), quantile=out.getQuantileExpected())
    )
    df.loc[df.index[:-1], "uncertainty"] = np.nan

    return df


def bayesian_simple(datacard, **kwargs):

    limit = combine(datacard, method="BayesianSimple", **kwargs).getLimit()[-1]

    return limit


def bayesian_toy_mc(datacard, **kwargs):

    limit = combine(datacard, method="BayesianToyMC", **kwargs).getLimit()[-1]

    return limit


def fit_diagnostics(datacard, **kwargs):

    import numpy as np
    import pandas as pd

    out = combine(datacard, method="FitDiagnostics", **kwargs)

    df = pd.DataFrame(
        data=dict(limit=out.getLimit(), uncertainty=out.getLimitError(), quantile=out.getQuantileExpected())
    )
    df.loc[df.index[:-1], "uncertainty"] = np.nan

    return df


def hybrid_new(datacard, **kwargs):

    res = combine(datacard, method="HybridNew", **kwargs)

    return res.getLimit()[-1], res.getLimitError()[-1]


def multi_dim_fit(datacard, algo=None, **kwargs):
     
    algos = {None : 0,
            "singles" : 1,
            "cross" : 2,
            "grid" : 3,
            "grid3x3" : 3, # note gridType_ should be G3x3 here
            "random" : 4,
            "contour2d" : 5,
            "stitch2d" : 6,
            "fixed" : 7,
            "impact" : 8,
            }

    if algo in algos:
        algo = algos[algo]
    else:
        raise ValueError("Unknown algorithm: " + algo)

    res = combine(datacard, method="MultiDimFit", algo=algo, **kwargs)

    return res.getR(), res.getDeltaNLL()


def significance(datacard, pvalue=False, **kwargs):
    import scipy.stats

    significance = combine(datacard, method="Significance", **kwargs).getLimit()[-1]

    if pvalue:
        return scipy.stats.norm.sf(significance)

    return significance
