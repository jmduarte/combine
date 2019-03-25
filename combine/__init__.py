import ROOT

from combine.DatacardParser import *
from combine.ModelTools import *
from combine.ShapeTools import *
from combine.PhysicsModel import *

from sys import modules


def text2workspace(
    fileName,
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
    nuisancesToExclude=[],
    nuisancesToRescale=[],
    nuisanceFunctions=[],
    nuisanceGroupFunctions=[],
    forceNonSimPdf=False,
    noCheckNorm=False,
    noJMax=False,
    allowNoSignal=False,
    allowNoBackground=False,
    optimizeExistingTemplates=True,
    optimizeBoundNuisances=True,
    optimizeTemplateBins=True,
    physModel="combine.PhysicsModel:defaultModel",
    physOpt=[],
    dumpCard=False,
    evaluateEdits=True,
):

    if fileName.endswith(".gz"):
        import gzip

        file = gzip.open(fileName, "rb")
        fileName = fileName[:-3]
    else:
        file = open(fileName, "r")

    ## Parse text file
    DC = parseCard(
        file, bin, noJMax, allowNoSignal, allowNoBackground, evaluateEdits, nuisancesToExclude, stat, verbose
    )

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
            fileName,
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
            fileName,
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
    libraries=[],
    compiled=False,
    keyword_value=[],
    text2workspace_kwargs={},
):
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
    )

    os.remove(tmp_file)

    return out


def significance(datacard, pvalue=False):
    import scipy.stats

    significance = combine(
        datacard, method="Significance", mass=125, toys=-1, expect_signal=1, significance=True
    ).getLimit()[-1]

    if pvalue:
        return scipy.stats.norm.sf(significance)

    return significance
