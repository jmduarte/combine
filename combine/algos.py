import ROOT
import numpy as np

from .DatacardParser import *
from .ModelTools import *
from .ShapeTools import *
from .PhysicsModel import *
from . import combine

from sys import exit
from optparse import OptionParser

import tempfile
import os


def binned_shape_significance(processes, signal_process, low, up, nbins, integrals=dict()):
    """Calculate the significance of a binned shape analysis without any systematics.

    The equivalent combine command would be:
        combine -n SignifExp -M Significance --significance -m 125 workspace.root -t -1 --expectSignal=1

    Args:
        processes (dict of key str and value array-like): A 1-d array of observables for different processes.
        signal_process (str): The  name of the process with represents the signal.
        low (float): Lower range limit for the histograms.
        up (float): Upper range limit for the histograms.
        nbin (int): Number of bins for the histograms.

    Returns:
        float: the significance value obtained by combine.
    """

    hist_file_name = os.path.join(tempfile.gettempdir(), next(tempfile._get_candidate_names()) + ".root")
    roof_file_name = os.path.join(tempfile.gettempdir(), next(tempfile._get_candidate_names()) + ".root")

    hists = dict()
    hists["data_obs"] = ROOT.TH1F("data_obs", "data_obs", nbins, low, up)

    weights = {k: float(v) / len(processes[k]) for k, v in integrals.items()}

    for k, v in processes.items():
        hists[k] = ROOT.TH1F(k, k, nbins, low, up)

    for k, v in processes.items():
        if k in weights:
            w = weights[k]
            for x in v:
                hists[k].Fill(x, w)
                hists["data_obs"].Fill(x, w)
        else:
            for x in v:
                hists[k].Fill(x)
                hists["data_obs"].Fill(x)

    integrals = {}

    myfile = ROOT.TFile(hist_file_name, "RECREATE")
    for k, v in hists.items():
        integrals[k] = v.Integral()
        v.Write()
    myfile.Close()

    del hists

    parser = OptionParser()
    addDatacardParserOptions(parser)
    options, args = parser.parse_args()
    options.bin = True  # make a binary workspace

    DC = Datacard()
    MB = None

    ############## Setup the datacard (must be filled in) ###########################

    DC.bins = ["bin1"]
    DC.obs = {"bin1": integrals["data_obs"]}
    DC.processes = list(processes.keys())
    DC.signals = [signal_process]
    DC.isSignal = {p: p == signal_process for p in DC.processes}
    DC.keyline = [("bin1", k, v) for k, v in DC.isSignal.items()]
    DC.exp = {"bin1": {k: v for k, v in integrals.items() if not k == "data_obs"}}
    DC.systs = []
    DC.shapeMap = {"*": {"*": [hist_file_name, "$PROCESS", "$PROCESS_$SYSTEMATIC"]}}
    DC.hasShapes = True
    DC.flatParamNuisances = {}
    DC.rateParams = {}
    DC.extArgs = {}
    DC.rateParamsOrder = set()
    DC.frozenNuisances = set()
    DC.systematicsShapeMap = {}
    DC.nuisanceEditLines = []
    DC.binParFlags = {}
    DC.groups = {}
    DC.discretes = []

    ###### User defined options #############################################

    options.out = roof_file_name  # Output workspace name
    options.fileName = "./"  # Path to input ROOT files
    options.verbose = 1  # Verbosity

    ##########################################################################

    if DC.hasShapes:
        MB = ShapeBuilder(DC, options)
    else:
        MB = CountingModelBuilder(DC, options)

    # Set physics models
    MB.setPhysics(defaultModel)
    MB.doModel()

    res = combine(
        roof_file_name, method="Significance", mass=125, toys=-1, expect_signal=1, significance=True, systematics=0
    )
    os.remove(hist_file_name)
    os.remove(roof_file_name)

    return res
