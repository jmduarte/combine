/**************************************
  Simple multiChannel significance & limit calculator
***************************************/
#include "combine/Combine.h"
#include <cstring>
#include <cerrno>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <unistd.h>
#include <errno.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TIterator.h>
#include <TLine.h>
#include <TMath.h>
#include <TString.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include <TTree.h>

#include <RooAbsData.h>
#include <RooAbsPdf.h>
#include <RooArgSet.h>
#include <RooCustomizer.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooMsgService.h>
#include <RooPlot.h>
#include <RooRandom.h>
#include <RooRealVar.h>
#include <RooUniform.h>
#include <RooGaussian.h>
#include <RooWorkspace.h>
#include <RooCategory.h>

#include <RooStats/HLFactory.h>
#include <RooStats/RooStatsUtils.h>
#include <RooStats/ModelConfig.h>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string.hpp>
#include <regex>

#include "combine/LimitAlgo.h"
#include "combine/utils.h"
#include "combine/CloseCoutSentry.h"
#include "combine/RooSimultaneousOpt.h"
#include "combine/ToyMCSamplerOpt.h"
#include "combine/AsimovUtils.h"
#include "combine/CascadeMinimizer.h"
#include "combine/ProfilingTools.h"

#include "combine/Logger.h"

using namespace RooStats;
using namespace RooFit;
using namespace std;

std::unique_ptr<LimitAlgo> algo;
std::unique_ptr<LimitAlgo> g_hintAlgo;

Float_t g_quantileExpected = -1.0;
TDirectory *g_outputFile = 0;
TDirectory *g_writeToysHere = 0;
TDirectory *g_readToysFromHere = 0;
int g_verbose = 1;
bool g_withSystematics = 1;
bool expectSignalSet_ = false;
bool g_doSignificance = 0;
bool expectSignalSet = 0;
bool g_lowerLimit = 0;
float g_confidenceLevel = 0.95;
bool g_bypassFrequentistFit = false;
bool g_fillTree = true;
TTree *Combine::tree_ = 0;
float g_mass = 120;

std::string g_setPhysicsModelParameterExpression = "";
std::string g_setPhysicsModelParameterRangeExpression = "";
std::string g_defineBackgroundOnlyModelParameterExpression = "";

std::string Combine::trackParametersNameString_ = "";
std::string Combine::textToWorkspaceString_ = "";

std::vector<std::pair<RooAbsReal *, float> > Combine::trackedParametersMap_;

Combine::Combine(float expectSignal, bool frequentistToys, bool toysNoSystematics)
    : statOptions_("Common statistics options"),
      ioOptions_("Common input-output options"),
      miscOptions_("Common miscellaneous options"),
      rMin_(std::numeric_limits<float>::quiet_NaN()),
      rMax_(std::numeric_limits<float>::quiet_NaN()),
      expectSignal_(expectSignal),
      toysFrequentist_(frequentistToys),
      toysNoSystematics_((!g_withSystematics) || toysNoSystematics) // if no systematics, also don't expect them for the toys
{
  namespace po = boost::program_options;
  statOptions_.add_options()("cl,C", po::value<float>(&g_confidenceLevel)->default_value(0.95), "Confidence Level")(
      "rMin", po::value<float>(&rMin_), "Override minimum value for signal strength (default is 0)")(
      "rMax", po::value<float>(&rMax_), "Override maximum value for signal strength (default is 20)")(
      "prior",
      po::value<std::string>(&prior_)->default_value("flat"),
      "Prior to use, for methods that require it and if it's not already in the input file: 'flat' (default), "
      "'1/sqrt(r)'")(
      "expectSignalMass",
      po::value<float>(&expectSignalMass_)->default_value(-99.),
      "If set to non-zero, generate *signal* toys instead of background ones, with the specified mass.")(
      "unbinned,U", "Generate unbinned datasets instead of binned ones (works only for extended pdfs)")(
      "generateBinnedWorkaround",
      "Make binned datasets generating unbinned ones and then binnning them. Workaround for a bug in RooFit.")(
      "setParameters",
      po::value<string>(&g_setPhysicsModelParameterExpression)->default_value(""),
      "Set the values of relevant physics model parameters. Give a comma separated list of parameter value "
      "assignments. Example: CV=1.0,CF=1.0")(
      "defineBackgroundOnlyModelParameters",
                               po::value<string>(&g_defineBackgroundOnlyModelParameterExpression)->default_value(""),
                               "If no background only (null) model is explicitly provided in physics model, one will "
                               "be defined as these values of the POIs (default is r=0)")(
      "redefineSignalPOIs",
      po::value<string>(&redefineSignalPOIs_)->default_value(""),
      "Redefines the POIs to be this comma-separated list of variables from the workspace.")(
      "freezeParameters",
      po::value<string>(&freezeNuisances_)->default_value(""),
      "Set as constant all these parameters.")("freezeNuisanceGroups",
                                               po::value<string>(&freezeNuisanceGroups_)->default_value(""),
                                               "Set as constant all these groups of nuisance parameters.");
  ioOptions_.add_options()("saveWorkspace", "Save workspace to output root file")(
      "workspaceName,w",
      po::value<std::string>(&workspaceName_)->default_value("w"),
      "Workspace name, when reading it from or writing it to a rootfile.")(
      "snapshotName",
      po::value<std::string>(&snapshotName_)->default_value(""),
      "Default snapshot name for pre-fit snapshot for reading or writing to workspace")(
      "modelConfigName",
      po::value<std::string>(&modelConfigName_)->default_value("ModelConfig"),
      "ModelConfig name, when reading it from or writing it to a rootfile.")(
      "modelConfigNameB",
      po::value<std::string>(&modelConfigNameB_)->default_value("%s_bonly"),
      "Name of the ModelConfig for b-only hypothesis.\n"
      "If not present, it will be made from the singal model taking zero signal strength.\n"
      "A '%s' in the name will be replaced with the modelConfigName.")(
      "overrideSnapshotMass", "Override MH loaded from a snapshot with the one passed on the command line")

      ("validateModel,V", "Perform some sanity checks on the model and abort if they fail.")(
          "saveToys", "Save results of toy MC in output file")(
          "floatAllNuisances",
          po::value<bool>(&floatAllNuisances_)->default_value(false),
          "Make all nuisance parameters floating")(
          "floatParameters",
          po::value<string>(&floatNuisances_)->default_value(""),
          "Set to floating these parameters (note freeze will take priority over float)")(
          "freezeAllGlobalObs",
          po::value<bool>(&freezeAllGlobalObs_)->default_value(true),
          "Make all global observables constant");
  miscOptions_.add_options()("optimizeSimPdf",
                             po::value<bool>(&optSimPdf_)->default_value(true),
                             "Turn on special optimizations of RooSimultaneous. On by default, you can turn it off if "
                             "it doesn't work for your workspace.")(
      "noMCbonly", po::value<bool>(&noMCbonly_)->default_value(false), "Don't create a background-only modelConfig")(
      "noDefaultPrior",
      po::value<bool>(&noDefaultPrior_)->default_value(false),
      "Don't create a default uniform prior")(
      "rebuildSimPdf",
      po::value<bool>(&rebuildSimPdf_)->default_value(false),
      "Rebuild simultaneous pdf from scratch to make sure constraints are correct (not needed in CMS workspaces)")(
      "compile", "Compile expressions instead of interpreting them")(
      "guessGenMode", "Guess if to generate binned or unbinned based on dataset")(
      "genBinnedChannels",
      po::value<std::string>(&genAsBinned_)->default_value(genAsBinned_),
      "Flag the given channels to be generated binned (irrespectively of how they were flagged at workspace creation)")(
      "genUnbinnedChannels",
      po::value<std::string>(&genAsUnbinned_)->default_value(genAsUnbinned_),
      "Flag the given channels to be generated unbinned (irrespectively of how they were flagged at workspace "
      "creation)")("text2workspace",
                   boost::program_options::value<std::string>(&textToWorkspaceString_)->default_value(""),
                   "Pass along options to text2workspace (default = none)")(
      "trackParameters",
      boost::program_options::value<std::string>(&trackParametersNameString_)->default_value(""),
      "Keep track of parameters in workspace, also accepts regexp with syntax 'rgx{<my regexp>}' (default = none)");
}

void Combine::applyOptions(std::string const &method, const boost::program_options::variables_map &vm) {
  if (g_withSystematics) {
    std::cout << ">>> including systematics" << std::endl;
  } else {
    std::cout << ">>> no systematics included" << std::endl;
  }
  unbinned_ = vm.count("unbinned");
  generateBinnedWorkaround_ = vm.count("generateBinnedWorkaround");
  if (unbinned_ && generateBinnedWorkaround_)
    throw std::logic_error("You can't set generateBinnedWorkaround and unbinned options at the same time");
  guessGenMode_ = vm.count("guessGenMode");
  compiledExpr_ = vm.count("compile");
  hintUsesStatOnly_ = vm.count("hintStatOnly");
  saveWorkspace_ = vm.count("saveWorkspace");
  if (modelConfigNameB_.find("%s") != std::string::npos) {
    char modelBName[1024];
    sprintf(modelBName, modelConfigNameB_.c_str(), modelConfigName_.c_str());
    modelConfigNameB_ = modelBName;
  }
  overrideSnapshotMass_ = vm.count("overrideSnapshotMass");
  saveToys_ = vm.count("saveToys");
  validateModel_ = vm.count("validateModel");
  if (!(vm["expectSignal"].defaulted()))
    expectSignalSet_ = true;
  else
    expectSignalSet_ = false;

  if (method == "MultiDimFit" || (method == "FitDiagnostics" && vm.count("justFit")) || method == "MarkovChainMC") {
    //CMSDAS new default,
    if (vm["noMCbonly"].defaulted())
      noMCbonly_ = 1;
    if (vm["noDefaultPrior"].defaulted())
      noDefaultPrior_ = 1;
  }
  if (!vm["prior"].defaulted())
    noDefaultPrior_ = 0;

  if (vm.count("keyword-value")) {
    modelPoints_ = vm["keyword-value"].as<std::vector<std::string> >();
  }
}

bool Combine::mklimit(RooWorkspace *workspace,
                      RooStats::ModelConfig *mc_s,
                      RooStats::ModelConfig *mc_b,
                      RooAbsData &data,
                      double &limit,
                      double &limitErr) {
  TStopwatch timer;

  bool ret = false;
  try {
    double hint = 0, hintErr = 0;
    bool hashint = false;
    if (g_hintAlgo) {
      if (hintUsesStatOnly_ && g_withSystematics) {
        g_withSystematics = false;
        hashint = g_hintAlgo->run(workspace, mc_s, mc_b, data, hint, hintErr, 0);
        g_withSystematics = true;
      } else {
        hashint = g_hintAlgo->run(workspace, mc_s, mc_b, data, hint, hintErr, 0);
      }
      workspace->loadSnapshot("clean");
    }
    limitErr = 0;  // start with 0, as some algorithms don't compute it
    ret = algo->run(workspace, mc_s, mc_b, data, limit, limitErr, (hashint ? &hint : 0));
  } catch (std::exception &ex) {
    std::cerr << "Caught exception " << ex.what() << std::endl;
    return false;
  }
  /* if ((ret == false) && (g_verbose > 3)) {
    std::cout << "Failed for method " << algo->name() << "\n";
    std::cout << "  --- DATA ---\n";
    utils::printRAD(&data);
    std::cout << "  --- MODEL ---\n";
    workspace->Print("V");
  } */
  timer.Stop();
  return ret;
}

void Combine::run(
    TString fileToLoad, const std::string &dataset, double &limit, double &limitErr, int &iToy, TTree *tree, int nToys) {
  ToCleanUp garbageCollect;  // use this to close and delete temporary files

  if (!boost::filesystem::exists(fileToLoad.Data()))
    throw std::invalid_argument(("File " + fileToLoad + " does not exist").Data());

  if (getenv("CMSSW_BASE")) {
    gSystem->AddIncludePath(TString::Format(" -I%s/src ", getenv("CMSSW_BASE")));
    if (g_verbose > 3)
      std::cout << "Adding " << getenv("CMSSW_BASE") << "/src to include path" << std::endl;
  }
  if (getenv("ROOFITSYS")) {
    gSystem->AddIncludePath(TString::Format(" -I%s/include ", getenv("ROOFITSYS")));
    if (g_verbose > 3)
      std::cout << "Adding " << getenv("ROOFITSYS") << "/include to include path" << std::endl;
  }

  if (g_verbose <= 2)
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  // Load the model, but going in a temporary directory to avoid polluting the current one with garbage from 'cexpr'
  RooWorkspace *workspace = 0;
  RooStats::ModelConfig *mc = 0, *mc_bonly = 0;
  std::unique_ptr<RooStats::HLFactory> hlf = nullptr;

  TFile *fIn = TFile::Open(fileToLoad);
  garbageCollect.tfile = fIn;  // request that we close this file when done

  workspace = dynamic_cast<RooWorkspace *>(fIn->Get(workspaceName_.c_str()));
  if (workspace == nullptr) {
    std::cerr << "Could not find workspace '" << workspaceName_ << "' in file " << fileToLoad << std::endl;
    fIn->ls();
    throw std::invalid_argument("Missing Workspace");
  }

  if (g_verbose > 3) {
    std::cout << "Input workspace '" << workspaceName_ << "': \n";
    workspace->Print("V");
  }
  RooRealVar *MH = workspace->var("MH");
  if (MH != 0) {
    if (g_verbose > 2)
      std::cerr << "Setting variable 'MH' in workspace to the higgs mass " << g_mass << std::endl;
    MH->setVal(g_mass);
  }
  mc = dynamic_cast<RooStats::ModelConfig *>(workspace->genobj(modelConfigName_.c_str()));
  mc_bonly = dynamic_cast<RooStats::ModelConfig *>(workspace->genobj(modelConfigNameB_.c_str()));

  if (mc == 0) {
    std::cerr << "Could not find ModelConfig '" << modelConfigName_ << "' in workspace '" << workspaceName_
              << "' in file " << fileToLoad << std::endl;
    throw std::invalid_argument("Missing ModelConfig");
  } else if (g_verbose > 3) {
    std::cout << "Workspace has a ModelConfig for signal called '" << modelConfigName_ << "', with contents:\n";
    mc->Print("V");
  }
  if (g_verbose > 3) {
    std::cout << "Input ModelConfig '" << modelConfigName_ << "': \n";
    mc->Print("V");
  }
  const RooArgSet *POI = mc->GetParametersOfInterest();
  if (POI == 0 || POI->getSize() == 0)
    throw std::invalid_argument("ModelConfig '" + modelConfigName_ + "' does not define parameters of interest.");
  if (POI->getSize() > 1)
    std::cerr << "ModelConfig '" << modelConfigName_
              << "' defines more than one parameter of interest. This is not supported in some statistical methods."
              << std::endl;
  if (mc->GetObservables() == 0)
    throw std::invalid_argument("ModelConfig '" + modelConfigName_ + "' does not define observables.");
  if (mc->GetPdf() == 0)
    throw std::invalid_argument("ModelConfig '" + modelConfigName_ + "' does not define a pdf.");
  if (rebuildSimPdf_ && typeid(*mc->GetPdf()) == typeid(RooSimultaneous)) {
    RooSimultaneous *newpdf =
        utils::rebuildSimPdf(*mc->GetObservables(), dynamic_cast<RooSimultaneous *>(mc->GetPdf()));
    workspace->import(*newpdf);
    mc->SetPdf(*newpdf);
  }
  if (optSimPdf_ && typeid(*mc->GetPdf()) == typeid(RooSimultaneous)) {
    RooSimultaneousOpt *optpdf = new RooSimultaneousOpt(static_cast<RooSimultaneous &>(*mc->GetPdf()),
                                                        TString(mc->GetPdf()->GetName()) + "_opt");
    workspace->import(*optpdf);
    mc->SetPdf(*optpdf);
  }
  if (mc_bonly == 0 && !noMCbonly_) {
    std::cerr << "Missing background ModelConfig '" << modelConfigNameB_ << "' in workspace '" << workspaceName_
              << "' in file " << fileToLoad << std::endl;
    RooCustomizer make_model_s(*mc->GetPdf(), "_model_bonly_");

    if (g_defineBackgroundOnlyModelParameterExpression != "") {
      std::cerr << "Will make one from the signal ModelConfig " << modelConfigName_ << " setting " << std::endl;
      vector<string> SetParameterExpressionList;
      boost::split(SetParameterExpressionList, g_defineBackgroundOnlyModelParameterExpression, boost::is_any_of(","));
      for (UInt_t p = 0; p < SetParameterExpressionList.size(); ++p) {
        vector<string> SetParameterExpression;
        boost::split(SetParameterExpression, SetParameterExpressionList[p], boost::is_any_of("="));
        if (SetParameterExpression.size() != 2) {
          std::cout << "Error parsing background model parameter expression : " << SetParameterExpressionList[p]
                    << endl;
        } else {
          std::string expr = SetParameterExpression[0];
          double expval = atof(SetParameterExpression[1].c_str());
          workspace->factory(Form("_%s_background_only_[%g]", expr.c_str(), expval));

          make_model_s.replaceArg(*POI->selectByName(expr.c_str())->first(),
                                  *workspace->var(Form("_%s_background_only_", expr.c_str())));
          std::cerr << "   " << expr << " to " << expval << std::endl;
        }
      }
    } else {
      std::cerr << "Will make one from the signal ModelConfig '" << modelConfigName_ << "' setting signal strenth '"
                << POI->first()->GetName() << "' to zero" << std::endl;
      workspace->factory("_zero_[0]");
      make_model_s.replaceArg(*POI->first(), *workspace->var("_zero_"));
    }

    RooAbsPdf *model_b = dynamic_cast<RooAbsPdf *>(make_model_s.build());
    model_b->SetName("_model_bonly_");
    workspace->import(*model_b);
    mc_bonly = new RooStats::ModelConfig(*mc);
    mc_bonly->SetPdf(*model_b);
  }

  // Specific settings should be executed before user specified ranges!
  RooRealVar *r = (RooRealVar *)POI->first();
  if (!isnan(rMin_))
    r->setMin(rMin_);
  if (!isnan(rMax_))
    r->setMax(rMax_);
  if (!isnan(rMin_) || !isnan(rMax_)) {
    r->setVal(0.5 * (r->getMin() + r->getMax()));
  }

  if (snapshotName_ != "") {
    bool loaded = workspace->loadSnapshot(snapshotName_.c_str());
    assert(loaded);
    if (MH) {
      //make sure mass value used is really the one from the loaded snapshot unless explicitly requested to override it
      if (overrideSnapshotMass_) {
        MH->setVal(g_mass);
      } else {
        g_mass = MH->getVal();
      }
    }
  }

  if (g_setPhysicsModelParameterRangeExpression != "") {
    utils::setModelParameterRanges(g_setPhysicsModelParameterRangeExpression, workspace->allVars());
  }
  //*********************************************
  //set physics model parameters    after loading the snapshot
  //*********************************************
  if (g_setPhysicsModelParameterExpression != "") {
    RooArgSet allParams(workspace->allVars());
    //if (workspace->genobj("discreteParams")) allParams.add(*(RooArgSet*)workspace->genobj("discreteParams"));
    allParams.add(workspace->allCats());
    utils::setModelParameters(g_setPhysicsModelParameterExpression, allParams);
    // also allow for "discrete" parameters to be set
    // Possible that MH value was re-set above, so make sure mass is set to the correct value and not over-ridden later.
    if (workspace->var("MH"))
      g_mass = workspace->var("MH")->getVal();
  }

  if (g_verbose <= 2)
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  const RooArgSet *observables = mc->GetObservables();       // not null
  POI = mc->GetParametersOfInterest();      // not null
  const RooArgSet *nuisances = mc->GetNuisanceParameters();  // note: may be null
  if (dynamic_cast<RooRealVar *>(POI->first()) == 0)
    throw std::invalid_argument("First parameter of interest is not a RooRealVar");

  if (nuisances && runtimedef::get("ADD_DISCRETE_FALLBACK")) {
    RooArgSet newNuis;
    std::string startswith = "u_CMS_Hgg_env_pdf_";
    TIterator *np = nuisances->createIterator();
    while (RooRealVar *arg = (RooRealVar *)np->Next()) {
      if (std::string(arg->GetName()).compare(0, startswith.size(), startswith)) {
        newNuis.add(*arg);
      } else {
        std::cout << "Removed nuisance from set: " << arg->GetName() << "\n";
      }
    }
    if (newNuis.getSize() < nuisances->getSize()) {
      mc->SetNuisanceParameters(newNuis);
      if (mc_bonly)
        mc_bonly->SetNuisanceParameters(newNuis);
      nuisances = mc->GetNuisanceParameters();
    }
  }

  if (dataset.find(":") != std::string::npos) {
    std::string filename, wspname, dname;
    switch (std::count(dataset.begin(), dataset.end(), ':')) {
      case 2:  // file:wsp:dataset
        filename = dataset.substr(0, dataset.find(":"));
        wspname = dataset.substr(dataset.find(":") + 1, dataset.rfind(":") - dataset.find(":") - 1);
        dname = dataset.substr(dataset.rfind(":") + 1, std::string::npos);
        if (g_verbose)
          std::cout << "Will read dataset '" << dname << "' from workspace '" << wspname << "' of file '" << filename
                    << "'" << std::endl;
        break;
      case 1:
        filename = dataset.substr(0, dataset.find(":"));
        dname = dataset.substr(dataset.find(":") + 1, std::string::npos);
        if (g_verbose)
          std::cout << "Will read dataset '" << dname << "' from file '" << filename << "'" << std::endl;
        break;
      default:
        throw std::invalid_argument("The dataset must be a name, or file:name or file:workspace:name");
    }
    if (filename == "" || dname == ":")
      throw std::invalid_argument("The dataset must be a name, or file:name or file:workspace:name");
    TDirectory *pwd = gDirectory;
    TFile *file = TFile::Open(filename.c_str());
    RooAbsData *data = 0;
    if (file == 0)
      throw std::invalid_argument(std::string("Cannot open input file: ") + filename);
    if (wspname.empty()) {
      data = (RooAbsData *)file->Get(dname.c_str());
      if (data == 0)
        throw std::invalid_argument(std::string("Cannot find a dataset named ") + dname + " in file " + filename);
    } else {
      RooWorkspace *win = (RooWorkspace *)file->Get(wspname.c_str());
      if (win == 0)
        throw std::invalid_argument(std::string("Cannot find a workspace named ") + wspname + " in file " + filename);
      data = (RooAbsData *)win->data(dname.c_str());
      if (data == 0)
        throw std::invalid_argument(std::string("Cannot find a dataset named ") + dname + " in file " + filename +
                                    ", workspace " + wspname);
    }
    workspace->import(*data, RooFit::Rename(dataset.c_str()));
    file->Close();
    pwd->cd();
  }
  if (workspace->data(dataset.c_str()) == 0) {
    TFile *fIn = TFile::Open(fileToLoad);
    garbageCollect.tfile = fIn;  // request that we close this file when done
    RooDataSet *data_obs = dynamic_cast<RooDataSet *>(fIn->Get(dataset.c_str()));
    if (data_obs) {
      data_obs->SetName(dataset.c_str());
      workspace->import(*data_obs);
    } else {
      std::cout << "Dataset " << dataset.c_str() << " not found." << std::endl;
    }
  }

  if (g_verbose < -1) {
    RooMsgService::instance().setStreamStatus(0, kFALSE);
    RooMsgService::instance().setStreamStatus(1, kFALSE);
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  }

  if (redefineSignalPOIs_ != "") {
    RooArgSet newPOIs(workspace->argSet(redefineSignalPOIs_.c_str()));
    TIterator *np = newPOIs.createIterator();
    while (RooRealVar *arg = (RooRealVar *)np->Next()) {
      RooRealVar *rrv = dynamic_cast<RooRealVar *>(arg);
      if (rrv == 0) {
        std::cerr << "combine/MultiDimFit: Parameter of interest " << arg->GetName()
                  << " which is not a RooRealVar will be ignored" << std::endl;
        continue;
      }
      arg->setConstant(0);
      // also set ignoreConstraint flag for constraint PDF
      if (workspace->pdf(Form("%s_Pdf", arg->GetName())))
        workspace->pdf(Form("%s_Pdf", arg->GetName()))->setAttribute("ignoreConstraint");
    }
    if (g_verbose > 0) {
      std::cout << "Redefining the POIs to be: ";
      newPOIs.Print("");
    }
    mc->SetParametersOfInterest(newPOIs);
    POI = mc->GetParametersOfInterest();
    if (nuisances) {
      RooArgSet newNuis(*nuisances);
      newNuis.remove(*POI);
      if (newNuis.getSize() < nuisances->getSize()) {
        mc->SetNuisanceParameters(newNuis);
        if (mc_bonly)
          mc_bonly->SetNuisanceParameters(newNuis);
        nuisances = mc->GetNuisanceParameters();
      }
    }
  }
  // Always reset the POIs to floating (post-fit workspaces can actually have them frozen in some cases, in any case they can be re-frozen in the next step
  TIterator *pois = POI->createIterator();
  while (RooRealVar *arg = (RooRealVar *)pois->Next()) {
    arg->setConstant(0);
  }

  if (floatNuisances_ != "") {
    RooArgSet toFloat((floatNuisances_ == "all") ? *nuisances : (workspace->argSet(floatNuisances_.c_str())));
    if (g_verbose > 0) {
      std::cout << "Set floating the following parameters: ";
      toFloat.Print("");
      Logger::instance().log(std::string(Form("Combine.cc: %d -- Set floating the following parameters: ", __LINE__)),
                             Logger::kLogLevelInfo,
                             __func__);
      std::unique_ptr<TIterator> iter(toFloat.createIterator());
      for (RooAbsArg *a = (RooAbsArg *)iter->Next(); a != 0; a = (RooAbsArg *)iter->Next()) {
        Logger::instance().log(
            std::string(Form("Combine.cc: %d  %s ", __LINE__, a->GetName())), Logger::kLogLevelInfo, __func__);
      }
    }
    utils::setAllConstant(toFloat, false);
  }

  if (freezeNuisances_ != "") {
    // expand regexps
    while (freezeNuisances_.find("rgx{") != std::string::npos) {
      size_t pos1 = freezeNuisances_.find("rgx{");
      size_t pos2 = freezeNuisances_.find("}", pos1);
      std::string prestr = freezeNuisances_.substr(0, pos1);
      std::string poststr = freezeNuisances_.substr(pos2 + 1, freezeNuisances_.size() - pos2);
      std::string reg_esp = freezeNuisances_.substr(pos1 + 4, pos2 - pos1 - 4);

      //std::cout<<"interpreting "<<reg_esp<<" as regex "<<std::endl;
      std::regex rgx(reg_esp, std::regex::ECMAScript);

      std::string matchingParams = "";
      std::unique_ptr<TIterator> iter(nuisances->createIterator());
      for (RooAbsArg *a = (RooAbsArg *)iter->Next(); a != 0; a = (RooAbsArg *)iter->Next()) {
        const std::string &target = a->GetName();
        std::smatch match;
        if (std::regex_match(target, match, rgx)) {
          matchingParams = matchingParams + target + ",";
        }
      }
      freezeNuisances_ = prestr + matchingParams + poststr;
      freezeNuisances_ = boost::replace_all_copy(freezeNuisances_, ",,", ",");
    }

    RooArgSet toFreeze((freezeNuisances_ == "all") ? *nuisances : (workspace->argSet(freezeNuisances_.c_str())));
    if (g_verbose > 0) {
      std::cout << "Freezing the following parameters: ";
      toFreeze.Print("");
      Logger::instance().log(std::string(Form("Combine.cc: %d -- Freezing the following parameters: ", __LINE__)),
                             Logger::kLogLevelInfo,
                             __func__);
      std::unique_ptr<TIterator> iter(toFreeze.createIterator());
      for (RooAbsArg *a = (RooAbsArg *)iter->Next(); a != 0; a = (RooAbsArg *)iter->Next()) {
        Logger::instance().log(
            std::string(Form("Combine.cc: %d  %s ", __LINE__, a->GetName())), Logger::kLogLevelInfo, __func__);
      }
    }
    utils::setAllConstant(toFreeze, true);
    if (nuisances) {
      RooArgSet newnuis(*nuisances);
      newnuis.remove(toFreeze, /*silent=*/true, /*byname=*/true);
      mc->SetNuisanceParameters(newnuis);
      if (mc_bonly)
        mc_bonly->SetNuisanceParameters(newnuis);
      nuisances = mc->GetNuisanceParameters();
    }
  }

  if (freezeNuisanceGroups_ != "") {
    std::vector<string> nuisanceGroups;
    boost::algorithm::split(nuisanceGroups, freezeNuisanceGroups_, boost::algorithm::is_any_of(","));
    for (std::vector<string>::iterator ng_it = nuisanceGroups.begin(); ng_it != nuisanceGroups.end(); ng_it++) {
      bool freeze_complement = false;
      if (boost::algorithm::starts_with((*ng_it), "^")) {
        freeze_complement = true;
        (*ng_it).erase(0, 1);
      }

      if (!workspace->set(Form("group_%s", (*ng_it).c_str()))) {
        std::cerr << "Unknown nuisance group: " << (*ng_it) << std::endl;
        throw std::invalid_argument("Unknown nuisance group name");
      }
      RooArgSet groupNuisances(*(workspace->set(Form("group_%s", (*ng_it).c_str()))));
      RooArgSet toFreeze;

      if (freeze_complement) {
        RooArgSet still_floating(*mc->GetNuisanceParameters());
        still_floating.remove(groupNuisances, true, true);
        toFreeze.add(still_floating);
      } else {
        toFreeze.add(groupNuisances);
      }

      if (g_verbose > 0) {
        std::cout << "Freezing the following nuisance parameters: ";
        toFreeze.Print("");
      }
      utils::setAllConstant(toFreeze, true);
      if (nuisances) {
        RooArgSet newnuis(*nuisances);
        newnuis.remove(toFreeze, /*silent=*/true, /*byname=*/true);
        mc->SetNuisanceParameters(newnuis);
        if (mc_bonly)
          mc_bonly->SetNuisanceParameters(newnuis);
        nuisances = mc->GetNuisanceParameters();
      }
    }
  }

  if (mc->GetPriorPdf() == 0 && !noDefaultPrior_) {
    if (prior_ == "flat") {
      RooAbsPdf *prior = new RooUniform("prior", "prior", *POI);
      workspace->import(*prior);
      mc->SetPriorPdf(*prior);
    } else if (prior_ == "1/sqrt(r)") {
      std::cout << "Will use prior 1/sqrt(" << POI->first()->GetName() << std::endl;
      TString priorExpr = TString::Format("EXPR::prior(\"1/sqrt(@0)\",%s)", POI->first()->GetName());
      workspace->factory(priorExpr.Data());
      mc->SetPriorPdf(*workspace->pdf("prior"));
    } else if (!prior_.empty() && workspace->pdf(prior_.c_str()) != 0) {
      std::cout << "Will use prior '" << prior_ << "' in from the input workspace" << std::endl;
      mc->SetPriorPdf(*workspace->pdf(prior_.c_str()));
    } else {
      std::cerr << "Unknown prior '" << prior_ << "'. It's not 'flat' '1/sqrt(r)' or the name of a pdf in the model.\n"
                << std::endl;
      throw std::invalid_argument("Bad prior");
    }
  }

  if (g_withSystematics && nuisances == 0) {
    std::cout << "The model has no constrained nuisance parameters. Please run the limit tool with no systematics "
                 "(option -S 0)."
              << std::endl;
    std::cout << "To make things easier, I will assume you have done it." << std::endl;
    if (g_verbose)
      Logger::instance().log(
          std::string(Form("Combine.cc: %d -- The signal model has no constrained nuisance parameters so I have "
                           "assumed you don't need a pdf for them. Please re-run with -S 0 to be sure!",
                           __LINE__)),
          Logger::kLogLevelInfo,
          __func__);
    g_withSystematics = false;
  } else if (!g_withSystematics && nuisances != 0) {
    std::cout << "Will set nuisance parameters to constants: ";
    utils::setAllConstant(*nuisances, true);
  }

  bool validModel = validateModel_ ? utils::checkModel(*mc, false) : true;
  if (validateModel_ && g_verbose)
    std::cout << "Sanity checks on the model: " << (validModel ? "OK" : "FAIL") << std::endl;

  // make sure these things are set consistently with what we expect
  if (floatAllNuisances_ && mc->GetNuisanceParameters() && g_withSystematics)
    utils::setAllConstant(*mc->GetNuisanceParameters(), false);
  if (freezeAllGlobalObs_ && mc->GetGlobalObservables())
    utils::setAllConstant(*mc->GetGlobalObservables(), true);
  if (floatAllNuisances_ && mc_bonly && mc_bonly->GetNuisanceParameters() && g_withSystematics)
    utils::setAllConstant(*mc_bonly->GetNuisanceParameters(), false);
  if (freezeAllGlobalObs_ && mc_bonly && mc_bonly->GetGlobalObservables())
    utils::setAllConstant(*mc_bonly->GetGlobalObservables(), true);

  // Setup the CascadeMinimizer with discrete nuisances
  addDiscreteNuisances(workspace);
  // and give him the regular nuisances too
  addNuisances(nuisances);
  addFloatingParameters(workspace->allVars());
  addPOI(POI);

  tree_ = tree;

  // Set up additional branches
  if (trackParametersNameString_ != "") {
    char tmp[10240];
    strlcpy(tmp, trackParametersNameString_.c_str(), 10240);
    char *token = strtok(tmp, ",");
    while (token) {
      if (boost::starts_with(token, "rgx{") && boost::ends_with(token, "}")) {
        std::string tokenstr(token);
        std::string reg_esp = tokenstr.substr(4, tokenstr.size() - 5);
        std::cout << "interpreting " << reg_esp << " as regex " << std::endl;
        std::regex rgx(reg_esp, std::regex::ECMAScript);

        RooArgSet allParams(workspace->allVars());
        std::unique_ptr<TIterator> iter(allParams.createIterator());
        for (RooAbsArg *a = (RooAbsArg *)iter->Next(); a != 0; a = (RooAbsArg *)iter->Next()) {
          RooAbsReal *tmp = dynamic_cast<RooAbsReal *>(a);
          const std::string &target = tmp->GetName();
          std::smatch match;
          if (std::regex_match(target, match, rgx)) {
            if (tmp->isConstant())
              continue;
            Combine::trackedParametersMap_.push_back(std::pair<RooAbsReal *, float>(tmp, tmp->getVal()));
          }
        }
        token = strtok(0, ",");
      } else {
        RooAbsReal *a = (RooAbsReal *)workspace->obj(token);
        if (a == 0)
          throw std::invalid_argument(std::string("Parameter ") + (token) + " not in model.");
        Combine::trackedParametersMap_.push_back(std::pair<RooAbsReal *, float>(a, a->getVal()));
        token = strtok(0, ",");
      }
    }
  }

  for (std::vector<std::pair<RooAbsReal *, float> >::iterator it = Combine::trackedParametersMap_.begin();
       it != Combine::trackedParametersMap_.end();
       it++) {
    const char *token = (it->first)->GetName();
    addBranch((std::string("trackedParam_") + token).c_str(),
              &(it->second),
              (std::string("trackedParam_") + token + std::string("/F")).c_str());
  }
  // Should have the PDF at this point, if not something is really odd?
  if (!(mc->GetPdf())) {
    std::cerr << " FATAL ERROR! PDF not found in ModelConfig. \n Try to build the workspace first with "
                 "text2workspace.py and run with the binary output."
              << std::endl;
    assert(0);
  }

  // Print list of channel masks
  RooSimultaneousOpt *simopt = dynamic_cast<RooSimultaneousOpt *>(mc->GetPdf());
  if (simopt && simopt->channelMasks().getSize() > 0) {
    int nChnMasks = simopt->channelMasks().getSize();
    int nActiveMasks = 0;
    for (int iMask = 0; iMask < nChnMasks; iMask++) {
      if (dynamic_cast<RooRealVar *>(simopt->channelMasks().at(iMask))->getVal() > 0)
        nActiveMasks++;
    }
    std::cout << ">>> " << nActiveMasks << " out of " << nChnMasks << " channels masked\n" << std::endl;
    if (g_verbose >= 2) {
      std::cout << ">>> Channel masks:\n";
      simopt->channelMasks().Print("v");
    }
  }

  bool isExtended = mc->GetPdf()->canBeExtended();
  MH = workspace->var("MH");
  RooAbsData *dobs = workspace->data(dataset.c_str());
  // Generate with signal model if r or other physics model parameters are defined
  RooAbsPdf *genPdf = (expectSignal_ > 0 || g_setPhysicsModelParameterExpression != "" || !mc_bonly)
                          ? mc->GetPdf()
                          : (mc_bonly ? mc_bonly->GetPdf() : 0);
  RooRealVar *weightVar_ = 0;  // will be needed for toy generation in some cases
  if (guessGenMode_ && genPdf && genPdf->InheritsFrom("RooSimultaneous") && (dobs != 0)) {
    utils::guessChannelMode(dynamic_cast<RooSimultaneous &>(*mc->GetPdf()), *dobs, g_verbose);
    if (mc_bonly)
      utils::guessChannelMode(dynamic_cast<RooSimultaneous &>(*mc_bonly->GetPdf()), *dobs, 0);
  }
  if (!genAsBinned_.empty() || !genAsUnbinned_.empty()) {
    RooSimultaneous *sim = dynamic_cast<RooSimultaneous *>(genPdf);
    if (!sim)
      throw std::invalid_argument(
          "Options genBinnedChannels and genUnbinnedChannels only work for RooSimultaneous pdfs");
    utils::setChannelGenModes(*sim, genAsBinned_, genAsUnbinned_, g_verbose);
  }
  // With the value of r we have to handle four cases:
  //   1) --expectSignal given and r appears in --setPhysicsModelParameters --> prefer expectSignal value but warn user
  //   2) --expectSignal given and r not in --setPhysicsModelParameters --> use --expectSignal value
  //   3) --expectSignal not given and r appears in --setPhysicsModelParameters --> use --setPhysicsModelParameters value
  //   4) --expectSignal not given and r not in --setPhysicsModelParameters --> use default --expectSignal value
  if (nToys != 0) {
    if (POI->find("r")) {
      // Was r also specified in --setPhysicsModelParameters?
      bool rInParamExp = false;
      if (g_setPhysicsModelParameterExpression != "") {
        vector<string> SetParameterExpressionList;
        boost::split(SetParameterExpressionList, g_setPhysicsModelParameterExpression, boost::is_any_of(","));
        for (UInt_t p = 0; p < SetParameterExpressionList.size(); ++p) {
          vector<string> SetParameterExpression;
          boost::split(SetParameterExpression, SetParameterExpressionList[p], boost::is_any_of("="));
          if (SetParameterExpression.size() == 2 && SetParameterExpression[0] == "r") {
            rInParamExp = true;
            break;
          }
        }
      }
      // Set the value of r in cases 1), 2) and 4)
      if (expectSignalSet_ || (!expectSignalSet_ && !rInParamExp && snapshotName_ == "")) {
        ((RooRealVar *)POI->find("r"))->setVal(expectSignal_);
      }
      if (expectSignalSet_ && rInParamExp) {
        std::cerr << "Warning: A value of r is specified in both the --setParameters "
                     "and --expectSignal options. The argument of --expectSignal will take "
                     "precedence\n";
      }
      if (MH && expectSignalMass_ > 0.) {
        MH->setVal(expectSignalMass_);
      }
    } else if (expectSignalSet_) {
      std::cerr << "Warning: option --expectSignal only applies to models with "
                   "the POI \"r\", use --setParameters to set the "
                   "values of the POIs for toy generation in this model\n";
    }
  }

  // Ok now we're ready to go lets save a "clean snapshot" for the current parameters state
  // workspace->allVars() misses the RooCategories, useful for some things - so need to include them. Set up a utils function for that
  workspace->saveSnapshot("clean", utils::returnAllVars(workspace));

  if (nToys <= 0) {  // observed or asimov
    workspace->saveSnapshot("toyGenSnapshot", utils::returnAllVars(workspace));
    iToy = nToys;
    if (iToy == -1) {
      if (g_readToysFromHere != 0) {
        dobs = dynamic_cast<RooAbsData *>(g_readToysFromHere->Get("toys/toy_asimov"));
        if (dobs == 0) {
          std::cerr << "Toy toy_asimov not found in " << g_readToysFromHere->GetName() << ". List follows:\n";
          g_readToysFromHere->ls();
          return;
        }
        if (toysFrequentist_ && mc->GetGlobalObservables()) {
          RooAbsCollection *snap =
              dynamic_cast<RooAbsCollection *>(g_readToysFromHere->Get("toys/toy_asimov_snapshot"));
          if (!snap) {
            std::cerr << "Snapshot of global observables toy_asimov_snapshot not found in "
                      << g_readToysFromHere->GetName() << ". List follows:\n";
            g_readToysFromHere->ls();
            return;
          }
          RooArgSet gobs(*mc->GetGlobalObservables());
          gobs.assignValueOnly(*snap);
          workspace->saveSnapshot("clean", utils::returnAllVars(workspace));
        }
      } else {
        if (genPdf == 0)
          throw std::invalid_argument(
              "You can't generate background-only toys if you have no background-only pdf in the workspace and you "
              "have set --noMCbonly");
        if (toysFrequentist_) {
          workspace->saveSnapshot("reallyClean", utils::returnAllVars(workspace));
          if (dobs == 0)
            throw std::invalid_argument("Frequentist Asimov datasets can't be generated without a real dataset to fit");
          RooArgSet gobsAsimov;
          utils::setAllConstant(*mc->GetParametersOfInterest(), true);  // Fix poi, before fit
          double poiVal = 0.;
          if (mc->GetParametersOfInterest()->getSize()) {
            poiVal = dynamic_cast<RooRealVar *>(mc->GetParametersOfInterest()->first())->getVal();
          }
          dobs = asimovutils::asimovDatasetWithFit(mc, *dobs, gobsAsimov, !g_bypassFrequentistFit, poiVal, g_verbose);
          if (mc->GetGlobalObservables()) {
            RooArgSet gobs(*mc->GetGlobalObservables());
            gobs = gobsAsimov;
          }
          utils::setAllConstant(*mc->GetParametersOfInterest(), false);
          workspace->saveSnapshot("clean", utils::returnAllVars(workspace));
        } else {
          toymcoptutils::SimPdfGenInfo newToyMC(*genPdf, *observables, !unbinned_);

          // print the values of the parameters used to generate the toy
          if (g_verbose > 2) {
            Logger::instance().log(
                std::string(Form("Combine.cc: %d -- Generate Asimov toy from parameter values ... ", __LINE__)),
                Logger::kLogLevelInfo,
                __func__);
            std::unique_ptr<TIterator> iter(genPdf->getParameters((const RooArgSet *)0)->createIterator());
            for (RooAbsArg *a = (RooAbsArg *)iter->Next(); a != 0; a = (RooAbsArg *)iter->Next()) {
              TString varstring = utils::printRooArgAsString(a);
              Logger::instance().log(std::string(Form("Combine.cc: %d -- %s", __LINE__, varstring.Data())),
                                     Logger::kLogLevelInfo,
                                     __func__);
            }
          }
          // Also save the current state of the tree here but specify the quantile as -2 (i.e not the default, something specific to the toys)
          if (saveToys_)
            commitPoint(false, -2);

          dobs = newToyMC.generateAsimov(weightVar_, g_verbose);  // as simple as that
        }
      }
    } else if (dobs == 0) {
      std::cerr << "No observed data '" << dataset << "' in the workspace. Cannot compute limit.\n" << std::endl;
      return;
    }
    if (saveToys_) {
      g_writeToysHere->WriteTObject(dobs, "toy_asimov");
      if (toysFrequentist_ && mc->GetGlobalObservables()) {
        RooAbsCollection *snap = mc->GetGlobalObservables()->snapshot();
        if (snap)
          g_writeToysHere->WriteTObject(snap, "toy_asimov_snapshot");
      }
    }
    std::cout << "Computing results starting from "
              << (iToy == 0 ? "observation (a-posteriori)" : "expected outcome (a-priori)") << std::endl;
    if (g_verbose)
      Logger::instance().log(
          std::string(Form("Combine.cc: %d -- Computing results starting from %s",
                           __LINE__,
                           (iToy == 0 ? "observation (a-posteriori)" : "expected outcome (a-priori)"))),
          Logger::kLogLevelInfo,
          __func__);
    if (MH)
      MH->setVal(g_mass);
    if (g_verbose > (isExtended ? 3 : 2))
      utils::printRAD(dobs);
    if (mklimit(workspace, mc, mc_bonly, *dobs, limit, limitErr))
      commitPoint(0, g_quantileExpected);  //tree->Fill();

    // Set the global flag to write output to the tree again since some Methods overwrite this to avoid the fill above.
    toggleGlobalFillTree(true);
  }

  std::vector<double> limitHistory;
  std::unique_ptr<RooAbsPdf> nuisancePdf;
  if (nToys > 0) {
    if (genPdf == 0)
      throw std::invalid_argument(
          "You can't generate background-only toys if you have no background-only pdf in the workspace and you have "
          "set --noMCbonly");
    toymcoptutils::SimPdfGenInfo newToyMC(*genPdf, *observables, !unbinned_);
    double expLimit = 0;
    unsigned int nLimits = 0;
    workspace->loadSnapshot("clean");
    RooDataSet *systDs = 0;
    RooArgSet allFloatingParameters = workspace->allVars();
    allFloatingParameters.remove(*mc->GetParametersOfInterest());
    int nFloatingNonPoiParameters = utils::countFloating(allFloatingParameters);
    if (nFloatingNonPoiParameters && !toysNoSystematics_ && (g_readToysFromHere == 0)) {
      if (nuisances == 0)
        throw std::logic_error(
            "Running with systematic variation in toys enabled, but I found floating parameters (which are not POIs) "
            "but no constrain terms have been defined in the datacard. If this is ok, re-run with -S 0");
      nuisancePdf.reset(utils::makeNuisancePdf(
          expectSignal_ || g_setPhysicsModelParameterExpression != "" || noMCbonly_ ? *mc : *mc_bonly));
      if (toysFrequentist_) {
        if (mc->GetGlobalObservables() == 0)
          throw std::logic_error("Cannot use toysFrequentist with no global observables");
        workspace->saveSnapshot("reallyClean", utils::returnAllVars(workspace));
        if (!g_bypassFrequentistFit) {
          utils::setAllConstant(*mc->GetParametersOfInterest(), true);
          if (dobs == 0)
            throw std::logic_error("Cannot use toysFrequentist with no input dataset");
          CloseCoutSentry sentry(g_verbose < 3);
          //genPdf->fitTo(*dobs, RooFit::Save(1), RooFit::Minimizer("Minuit2","minimize"), RooFit::Strategy(0), RooFit::Hesse(0), RooFit::Constrain(*(expectSignal_ ?mc:mc_bonly)->GetNuisanceParameters()));
          std::unique_ptr<RooAbsReal> nll(genPdf->createNLL(
              *dobs,
              RooFit::Constrain(
                  *(expectSignal_ || g_setPhysicsModelParameterExpression != "" || noMCbonly_ ? mc : mc_bonly)
                       ->GetNuisanceParameters()),
              RooFit::Extended(genPdf->canBeExtended())));
          CascadeMinimizer minim(*nll, CascadeMinimizer::Constrained);
          minim.setStrategy(1);
          minim.minimize();
          utils::setAllConstant(*mc->GetParametersOfInterest(), false);
          workspace->saveSnapshot("clean", utils::returnAllVars(workspace));
        }
        if (nuisancePdf.get())
          systDs = nuisancePdf->generate(*mc->GetGlobalObservables(), nToys);
      } else {
        if (nuisancePdf.get())
          systDs = nuisancePdf->generate(*nuisances, nToys);
      }
    }
    std::unique_ptr<RooArgSet> vars(genPdf->getVariables());
    algo->setNToys(nToys);

    for (iToy = 1; iToy <= nToys; ++iToy) {
      algo->setToyNumber(iToy - 1);
      RooAbsData *absdata_toy = 0;
      if (g_readToysFromHere == 0) {
        workspace->loadSnapshot("clean");
        if (g_verbose > 3)
          utils::printPdf(genPdf);
        if (g_withSystematics && !toysNoSystematics_) {
          if (systDs) {
            if (systDs->numEntries() >= iToy)
              *vars = *systDs->get(iToy - 1);
          }
          if (toysFrequentist_)
            workspace->saveSnapshot("clean", utils::returnAllVars(workspace));
          if (g_verbose > 3)
            utils::printPdf(genPdf);
        }
        /* No longer need to set this because "clean" state is already set correctly even without toysFrequentist
        if (POI->find("r")) {
          if (expectSignal_) ((RooRealVar*)POI->find("r"))->setVal(expectSignal_);
        }
        */
        std::cout << "Generate toy " << iToy << "/" << nToys << std::endl;
        if (g_verbose > 2) {
          Logger::instance().log(
              std::string(
                  Form("Combine.cc: %d -- Generating toy %d/%d, from parameter values ... ", __LINE__, iToy, nToys)),
              Logger::kLogLevelInfo,
              __func__);
          std::unique_ptr<TIterator> iter(genPdf->getParameters((const RooArgSet *)0)->createIterator());
          for (RooAbsArg *a = (RooAbsArg *)iter->Next(); a != 0; a = (RooAbsArg *)iter->Next()) {
            TString varstring = utils::printRooArgAsString(a);
            Logger::instance().log(
                std::string(Form("Combine.cc: %d -- %s", __LINE__, varstring.Data())), Logger::kLogLevelInfo, __func__);
          }
        }

        // Also save the current state of the tree here but specify the quantile as -2 (i.e not the default, something specific to the toys)
        if (saveToys_)
          commitPoint(false, -2);
        if (isExtended) {
          absdata_toy = newToyMC.generate(weightVar_);  // as simple as that
        } else {
          RooDataSet *data_toy = genPdf->generate(*observables, 1);
          absdata_toy = data_toy;
        }
      } else {
        workspace->loadSnapshot(
            "clean");  // (*) this is needed in case running over toys+fits, to avoid starting from previous fit
        //-- constraints are set to toy values if frequentist (below) or set back to 0 (unecessarily) here.
        absdata_toy = dynamic_cast<RooAbsData *>(g_readToysFromHere->Get(TString::Format("toys/toy_%d", iToy)));
        if (absdata_toy == 0) {
          std::cerr << "Toy toy_" << iToy << " not found in " << g_readToysFromHere->GetName() << ". List follows:\n";
          g_readToysFromHere->ls();
          return;
        }
        if (toysFrequentist_ && mc->GetGlobalObservables()) {
          RooAbsCollection *snap =
              dynamic_cast<RooAbsCollection *>(g_readToysFromHere->Get(TString::Format("toys/toy_%d_snapshot", iToy)));
          if (!snap) {
            std::cerr << "Snapshot of global observables toy_" << iToy << "_snapshot not found in "
                      << g_readToysFromHere->GetName() << ". List follows:\n";
            g_readToysFromHere->ls();
            return;
          }
          vars->assignValueOnly(*snap);
          // note, we save over the "clean" values also for the parameters, so we've made sure they are the same as they were in (*)
          workspace->saveSnapshot("clean", utils::returnAllVars(workspace));
        }
      }
      if (g_verbose > (isExtended ? 3 : 2))
        utils::printRAD(absdata_toy);
      if (!toysFrequentist_)
        workspace->saveSnapshot("toyGenSnapshot", utils::returnAllVars(workspace));
      workspace->loadSnapshot("clean");
      if (toysFrequentist_)
        workspace->saveSnapshot("toyGenSnapshot", utils::returnAllVars(workspace));
      //if (g_verbose > 1) utils::printPdf(workspace, "model_b");
      if (mklimit(workspace, mc, mc_bonly, *absdata_toy, limit, limitErr)) {
        commitPoint(0, g_quantileExpected);  //tree->Fill();
        ++nLimits;
        expLimit += limit;
        limitHistory.push_back(limit);
      }
      // Set the global flag to write output to the tree again since some Methods overwrite this to avoid the fill above.
      toggleGlobalFillTree(true);

      if (saveToys_) {
        g_writeToysHere->WriteTObject(absdata_toy, TString::Format("toy_%d", iToy));
        if (toysFrequentist_ && mc->GetGlobalObservables()) {
          RooAbsCollection *snap = mc->GetGlobalObservables()->snapshot();
          g_writeToysHere->WriteTObject(snap, TString::Format("toy_%d_snapshot", iToy));
        }
      }
      delete absdata_toy;
    }
    if (weightVar_)
      delete weightVar_;
    expLimit /= nLimits;
    double rms = 0;
    for (std::vector<double>::const_iterator itl = limitHistory.begin(); itl != limitHistory.end(); ++itl) {
      rms += pow(*itl - expLimit, 2);
    }
    if (nLimits > 1) {
      rms = sqrt(rms / (nLimits - 1) / nLimits);
      cout << "mean   expected limit: r < " << expLimit << " +/- " << rms << " @ " << g_confidenceLevel * 100 << "%CL ("
           << nLimits << " toyMC)" << endl;
    } else {
      cout << "mean   expected limit: r < " << expLimit << " @ " << g_confidenceLevel * 100 << "%CL (" << nLimits
           << " toyMC)" << endl;
    }
    sort(limitHistory.begin(), limitHistory.end());
    if (nLimits > 0) {
      double medianLimit = (nLimits % 2 == 0 ? 0.5 * (limitHistory[nLimits / 2 - 1] + limitHistory[nLimits / 2])
                                             : limitHistory[nLimits / 2]);
      cout << "median expected limit: r < " << medianLimit << " @ " << g_confidenceLevel * 100 << "%CL (" << nLimits
           << " toyMC)" << endl;
      double hi68 = limitHistory[min<int>(nLimits - 1, ceil(0.84 * nLimits))];
      double lo68 = limitHistory[min<int>(nLimits - 1, floor(0.16 * nLimits))];
      double hi95 = limitHistory[min<int>(nLimits - 1, ceil(0.975 * nLimits))];
      double lo95 = limitHistory[min<int>(nLimits - 1, floor(0.025 * nLimits))];
      cout << "   68% expected band : " << lo68 << " < r < " << hi68 << endl;
      cout << "   95% expected band : " << lo95 << " < r < " << hi95 << endl;
    }
  }

  if (saveWorkspace_) {
    workspace->SetName(workspaceName_.c_str());
    workspace->loadSnapshot("clean");
    g_outputFile->WriteTObject(workspace, workspaceName_.c_str());
  }
}

void Combine::toggleGlobalFillTree(bool flag) { g_fillTree = flag; }

void Combine::commitPoint(bool expected, float quantile) {
  Float_t saveQuantile = g_quantileExpected;
  g_quantileExpected = quantile;

  for (std::vector<std::pair<RooAbsReal *, float> >::iterator it = Combine::trackedParametersMap_.begin();
       it != Combine::trackedParametersMap_.end();
       it++) {
    it->second = (it->first)->getVal();
  }

  if (g_fillTree)
    tree_->Fill();
  g_quantileExpected = saveQuantile;
}

void Combine::addBranch(const char *name, void *address, const char *leaflist) {
  tree_->Branch(name, address, leaflist);
}
void Combine::addPOI(const RooArgSet *poi) {
  // RooArgSet *nuisances = (RooArgSet*) workspace->set("nuisances");
  CascadeMinimizerGlobalConfigs::O().parametersOfInterest = RooArgList();
  if (poi != 0) {
    TIterator *pp = poi->createIterator();
    while (RooAbsArg *arg = (RooAbsArg *)pp->Next())
      (CascadeMinimizerGlobalConfigs::O().parametersOfInterest).add(*arg);
  }
}

void Combine::addNuisances(const RooArgSet *nuisances) {
  // RooArgSet *nuisances = (RooArgSet*) workspace->set("nuisances");
  CascadeMinimizerGlobalConfigs::O().nuisanceParameters = RooArgList();
  if (nuisances != 0) {
    TIterator *np = nuisances->createIterator();
    while (RooAbsArg *arg = (RooAbsArg *)np->Next())
      (CascadeMinimizerGlobalConfigs::O().nuisanceParameters).add(*arg);
  }
}
void Combine::addFloatingParameters(const RooArgSet &parameters) {
  CascadeMinimizerGlobalConfigs::O().allFloatingParameters = RooArgList();
  //if (parameters != 0) {
  TIterator *np = parameters.createIterator();
  while (RooAbsArg *arg = (RooAbsArg *)np->Next()) {
    if (!arg->isConstant())
      (CascadeMinimizerGlobalConfigs::O().allFloatingParameters).add(*arg);
  }
  //}
}
void Combine::addDiscreteNuisances(RooWorkspace *workspace) {
  RooArgSet *discreteParameters = (RooArgSet *)workspace->genobj("discreteParams");

  CascadeMinimizerGlobalConfigs::O().pdfCategories = RooArgList();
  CascadeMinimizerGlobalConfigs::O().allRooMultiPdfParams = RooArgList();

  if (discreteParameters != 0) {
    TIterator *dp = discreteParameters->createIterator();
    while (RooAbsArg *arg = (RooAbsArg *)dp->Next()) {
      RooCategory *cat = dynamic_cast<RooCategory *>(arg);
      if (cat && (!cat->isConstant() || runtimedef::get("ADD_DISCRETE_FALLBACK"))) {
        if (g_verbose) {
          std::cout << "Adding discrete " << cat->GetName() << "\n";
          if (g_verbose)
            Logger::instance().log(std::string(Form("Combine.cc: %d -- Adding discrete %s ", __LINE__, cat->GetName())),
                                   Logger::kLogLevelInfo,
                                   __func__);
        }
        (CascadeMinimizerGlobalConfigs::O().pdfCategories).add(*arg);
      }
    }
  }
  // Run through all of the categories in the workspace and look for "pdfindex" -> fall back option
  else if (runtimedef::get("ADD_DISCRETE_FALLBACK")) {
    RooArgSet discreteParameters_C = workspace->allCats();
    TIterator *dp = discreteParameters_C.createIterator();
    while (RooAbsArg *arg = (RooAbsArg *)dp->Next()) {
      RooCategory *cat = dynamic_cast<RooCategory *>(arg);
      if (!(std::string(cat->GetName()).find("pdfindex") != std::string::npos))
        continue;
      if (cat /* && !cat->isConstant()*/) {
        if (g_verbose) {
          std::cout << "Adding discrete " << cat->GetName() << "\n";
          if (g_verbose)
            Logger::instance().log(std::string(Form("Combine.cc: %d -- Adding discrete %s ", __LINE__, cat->GetName())),
                                   Logger::kLogLevelInfo,
                                   __func__);
        }
        (CascadeMinimizerGlobalConfigs::O().pdfCategories).add(*arg);
      }
    }
  }
  // Now lets go through the list of parameters which are associated to this discrete nuisance
  RooArgSet clients;
  utils::getClients(CascadeMinimizerGlobalConfigs::O().pdfCategories, (workspace->allPdfs()), clients);
  TIterator *it = clients.createIterator();
  while (RooAbsArg *arg = (RooAbsArg *)it->Next()) {
    RooAbsPdf *pdf = dynamic_cast<RooAbsPdf *>(arg);
    RooArgSet *pdfPars = pdf->getParameters((const RooArgSet *)0);
    std::unique_ptr<TIterator> iter_v(pdfPars->createIterator());
    for (RooAbsArg *a = (RooAbsArg *)iter_v->Next(); a != 0; a = (RooAbsArg *)iter_v->Next()) {
      RooRealVar *v = dynamic_cast<RooRealVar *>(a);
      if (!v)
        continue;
      if (!(v->isConstant()))
        (CascadeMinimizerGlobalConfigs::O().allRooMultiPdfParams).add(*v);
    }
  }
}
