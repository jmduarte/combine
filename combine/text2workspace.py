#!/usr/bin/env python
from __future__ import absolute_import
import re
from sys import argv, stdout, stderr, exit, modules

# import ROOT with a fix to get batch mode (http://root.cern.ch/phpBB3/viewtopic.php?t=3198)
argv.append("-b-")
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
argv.remove("-b-")

from combine.DatacardParser import *
from combine.ModelTools import *
from combine.ShapeTools import *
from combine.PhysicsModel import *


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


text2workspace("simple-shapes-TH1.txt")
