#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include <string>

#include <TString.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <RooRandom.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <boost/program_options.hpp>
#include <map>

#include "combine/Combine.h"
#include "combine/Significance.h"
#include "combine/HybridNew.h"
#include "combine/BayesianFlatPrior.h"
#include "combine/BayesianToyMC.h"
#include "combine/MarkovChainMC.h"
#include "combine/FeldmanCousins.h"
#include "combine/FitDiagnostics.h"
#include "combine/AsymptoticLimits.h"
#include "combine/GoodnessOfFit.h"
#include "combine/ChannelCompatibilityCheck.h"
#include "combine/MultiDimFit.h"
#include "combine/CascadeMinimizer.h"
#include "combine/ProfilingTools.h"
#include "combine/GenerateOnly.h"
#include "combine/Logger.h"

template <class T> struct ValueWithError {
  T value;
  T error;
};

double combine(std::string const &datacard,
               std::string const &name,
               int argc,
               char **argv,
               std::string const &dataset,
               std::string const &whichMethod,
               std::string const &whichHintMethod,
               std::string const &toysFile,
               int runToys,
               float expectSignal,
               int seed,
               bool perfCounters) {
  using namespace std;
  using namespace boost;
  namespace po = boost::program_options;

  vector<string> runtimeDefines;
  vector<string> modelPoints;
  vector<string> modelParamNameVector_;
  vector<string> modelParamValVector_;

  Combine combiner(expectSignal);

  map<string, LimitAlgo *> methods;
  algo = new Significance();
  methods.insert(make_pair(algo->name(), algo));
  algo = new BayesianFlatPrior();
  methods.insert(make_pair(algo->name(), algo));
  algo = new BayesianToyMC();
  methods.insert(make_pair(algo->name(), algo));
  algo = new MarkovChainMC();
  methods.insert(make_pair(algo->name(), algo));
  algo = new HybridNew();
  methods.insert(make_pair(algo->name(), algo));
  algo = new FeldmanCousins();
  methods.insert(make_pair(algo->name(), algo));
  algo = new FitDiagnostics();
  methods.insert(make_pair(algo->name(), algo));
  algo = new AsymptoticLimits();
  methods.insert(make_pair(algo->name(), algo));
  algo = new GoodnessOfFit();
  methods.insert(make_pair(algo->name(), algo));
  algo = new ChannelCompatibilityCheck();
  methods.insert(make_pair(algo->name(), algo));
  algo = new MultiDimFit();
  methods.insert(make_pair(algo->name(), algo));
  algo = new GenerateOnly();
  methods.insert(make_pair(algo->name(), algo));

  CascadeMinimizer::initOptions();

  string methodsDesc("Method to extract upper limit. Supported methods are: ");
  for (map<string, LimitAlgo *>::const_iterator i = methods.begin(); i != methods.end(); ++i) {
    if (i != methods.begin())
      methodsDesc += ", ";
    methodsDesc += i->first;
  }

  po::options_description desc("Main options");
  combiner.miscOptions().add_options()(
      "keyword-value",
      po::value<vector<string> >(&modelPoints),
      "Set keyword values with 'WORD=VALUE', will replace $WORD with VALUE in datacards. Filename will also be "
      "extended with 'WORDVALUE'. Can specify multiple times")(
      "X-rtd",
      po::value<vector<string> >(&runtimeDefines),
      "Define some constants to be used at runtime (for debugging purposes). The syntax is --X-rtd identifier[=value], "
      "where value is an integer and defaults to 1. Can specify multiple times")(
      "X-fpeMask", po::value<int>(), "Set FPE mask: 1=NaN, 2=Div0, 4=Overfl, 8=Underf, 16=Inexact; 7=default");
  desc.add(combiner.statOptions());
  desc.add(combiner.ioOptions());
  desc.add(CascadeMinimizer::options());
  desc.add(combiner.miscOptions());
  po::positional_options_description p;
  po::variables_map vm, vm0;

  // parse the first time, using only common options and allow unregistered options
  po::store(po::command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm0);
  po::notify(vm0);

  // now search for algo, and add option
  map<string, LimitAlgo *>::const_iterator it_algo = methods.find(whichMethod);
  desc.add(it_algo->second->options());

  // parse the first time, now include options of the algo but not unregistered ones
  po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
  po::notify(vm);

  combiner.applyOptions(whichMethod, vm);
  CascadeMinimizer::applyOptions(vm);

  algo = it_algo->second;
  algo->applyOptions(vm);
  cout << ">>> method used is " << whichMethod << endl;

  if (!whichHintMethod.empty()) {
    g_hintAlgo = methods.at(whichHintMethod);
    g_hintAlgo->applyDefaultOptions();
    cout << ">>> method used to hint where the upper limit is " << whichHintMethod << endl;
  }

  if (seed == -1) {
    if (g_verbose > 0)
      std::cout << ">>> Using OpenSSL to get a really random seed " << std::endl;
    FILE *rpipe = popen("openssl rand 8", "r");
    std::cout << ">>> Used OpenSSL to get a really random seed " << seed << std::endl;
  } else {
    std::cout << ">>> random number generator seed is " << seed << std::endl;
  }
  RooRandom::randomGenerator()->SetSeed(seed);

  TString massName = TString::Format("mH%g.", g_mass);
  TString toyName = "";
  if (runToys > 0 || seed != 123456 || vm.count("saveToys"))
    toyName = TString::Format("%d.", seed);
  if (vm.count("expectedFromGrid") && !vm["expectedFromGrid"].defaulted())
    toyName += TString::Format("quant%.3f.", vm["expectedFromGrid"].as<float>());
  if (vm.count("expected") && !vm["expected"].defaulted())
    toyName += TString::Format("quant%.3f.", vm["expected"].as<float>());

  if (vm.count("keyword-value")) {
    for (vector<string>::const_iterator rmp = modelPoints.begin(), endrmp = modelPoints.end(); rmp != endrmp; ++rmp) {
      std::string::size_type idx = rmp->find('=');
      if (idx == std::string::npos) {
        cerr << "No value found for keyword :\n\t" << *rmp << " use --keyword-value WORD=VALUE " << std::endl;
      } else {
        std::string name = rmp->substr(0, idx);
        std::string svalue = rmp->substr(idx + 1);
        if (g_verbose > 0)
          std::cout << "Setting keyword " << name << " to " << svalue << std::endl;
        modelParamNameVector_.push_back(name);
        modelParamValVector_.push_back(svalue);
        massName += TString(name.c_str()) + TString(svalue.c_str()) + ".";
      }
    }
  }

  TString fileName = "higgsCombine" + name + "." + whichMethod + "." + massName + toyName + "root";
  TFile *test = new TFile(fileName, "RECREATE");
  g_outputFile = test;
  TTree *t = new TTree("limit", "limit");
  int syst, iToy, iChannel;
  ValueWithError<double> limit;
  t->Branch("limit", &limit.value, "limit/D");
  t->Branch("limitErr", &limit.error, "limitErr/D");
  t->Branch("syst", &syst, "syst/I");
  t->Branch("iToy", &iToy, "iToy/I");
  t->Branch("iChannel", &iChannel, "iChannel/I");
  t->Branch("quantileExpected", &g_quantileExpected, "quantileExpected/F");
  for (unsigned int mpi = 0; mpi < modelParamNameVector_.size(); ++mpi) {
    std::string name = modelParamNameVector_[mpi];
    t->Branch(Form("%s", name.c_str()), &modelParamValVector_[mpi]);
  }

  g_writeToysHere = test->mkdir("toys", "toys");
  if (toysFile != "")
    g_readToysFromHere = TFile::Open(toysFile.c_str());

  syst = g_withSystematics;
  iChannel = 0;

  if (vm.count("igpMem"))
    setupIgProfDumpHook();

  if (vm.count("X-fpeMask"))
    gSystem->SetFPEMask(vm["X-fpeMask"].as<int>());

  // CMSDAS Defaults (you can turn off with --X-rtd <name>=0
  runtimedef::set("OPTIMIZE_BOUNDS", 1);
  runtimedef::set("ADDNLL_RECURSIVE", 1);
  runtimedef::set("ADDNLL_GAUSSNLL", 1);
  runtimedef::set("ADDNLL_HISTNLL", 1);
  runtimedef::set("ADDNLL_CBNLL", 1);
  runtimedef::set("TMCSO_AdaptivePseudoAsimov", 1);
  // Optimization for bare RooFit likelihoods (--optimizeSimPdf=0)
  runtimedef::set("MINIMIZER_optimizeConst", 2);
  runtimedef::set("MINIMIZER_rooFitOffset", 1);
  // Optimization for ATLAS HistFactory likelihoods
  runtimedef::set("ADDNLL_ROOREALSUM_FACTOR", 1);
  runtimedef::set("ADDNLL_ROOREALSUM_NONORM", 1);
  runtimedef::set("ADDNLL_ROOREALSUM_BASICINT", 1);
  runtimedef::set("ADDNLL_ROOREALSUM_KEEPZEROS", 1);
  runtimedef::set("ADDNLL_PRODNLL", 1);
  runtimedef::set("ADDNLL_HFNLL", 1);
  runtimedef::set("ADDNLL_HISTFUNCNLL", 1);
  runtimedef::set("ADDNLL_ROOREALSUM_CHEAPPROD", 1);

  for (vector<string>::const_iterator rtdp = runtimeDefines.begin(), endrtdp = runtimeDefines.end(); rtdp != endrtdp;
       ++rtdp) {
    std::string::size_type idx = rtdp->find('=');
    if (idx == std::string::npos) {
      runtimedef::set(*rtdp, 1);
      if (g_verbose > 0)
        std::cout << "Turning on runtime-define " << *rtdp << std::endl;
    } else {
      std::string name = rtdp->substr(0, idx);
      std::string svalue = rtdp->substr(idx + 1);
      int ivalue = atoi(svalue.c_str());
      if (g_verbose > 0)
        std::cout << "Setting runtime-define " << name << " to " << ivalue << std::endl;
      runtimedef::set(name, ivalue);
    }
  }

  try {
    combiner.run(datacard, dataset, limit.value, limit.error, iToy, t, runToys);
    if (g_verbose > 0)
      Logger::instance().printLog();
  } catch (std::exception &ex) {
    cerr << "Error when running the combination:\n\t" << ex.what() << std::endl;
    test->Close();
  }

  test->WriteTObject(t);
  test->Close();

  for (auto i = methods.begin(); i != methods.end(); ++i)
    delete i->second;

  if (perfCounters)
    PerfCounter::printAll();

  return limit.value;
}

double _combine(std::string const &datacard,
                std::string const &name,
                std::string const &dataset,
                std::string const &method,
                std::string const &hintMethod,
                std::string const &toysFile,
                int verbose,
                bool significance,
                float mass,
                bool withSystematics,
                int toys,
                float expectSignal,
                int seed,
                bool lowerLimit,
                bool bypassFrequentistFit,
                bool perfCounters) {
  std::vector<char *> args;

  args.push_back(const_cast<char *>("combine"));

  g_verbose = verbose;
  g_mass = mass;
  g_withSystematics = withSystematics;
  g_doSignificance = significance;
  g_lowerLimit = lowerLimit;
  g_bypassFrequentistFit = bypassFrequentistFit;

  return combine(datacard,
                 name,
                 args.size(),
                 &args[0],
                 dataset,
                 method,
                 hintMethod,
                 toysFile,
                 toys,
                 expectSignal,
                 seed,
                 perfCounters);
}

PYBIND11_MODULE(_combine, m) {
  m.doc() = "CMS Higgs Combination toolkit.";

  m.def("_combine", &_combine, "CMS Higgs Combination toolkit.");
}
