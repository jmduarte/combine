#include "combine/MultiDimFit.h"
#include <stdexcept>
#include <cmath>

#include "TMath.h"
#include "TFile.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooRandom.h"
#include "RooAbsData.h"
#include "RooCategory.h"
#include "RooFitResult.h"
//#include "RooMinimizerOpt.h"
#include "RooMinimizer.h"
#include <RooStats/ModelConfig.h>
#include "combine/Combine.h"
#include "combine/CascadeMinimizer.h"
#include "combine/CloseCoutSentry.h"
#include "combine/utils.h"
#include "combine/RobustHesse.h"

#include <Math/Minimizer.h>
#include <Math/MinimizerOptions.h>
#include <Math/QuantFuncMathCore.h>
#include <Math/ProbFunc.h>

using namespace RooStats;

std::string MultiDimFit::name_ = "";
MultiDimFit::Algo MultiDimFit::g_algo = MultiDimFit::Algo::None;
MultiDimFit::GridType MultiDimFit::gridType_ = G1x1;
std::vector<std::string> MultiDimFit::poi_;
std::vector<RooRealVar *> MultiDimFit::poiVars_;
std::vector<float> MultiDimFit::poiVals_;
RooArgList MultiDimFit::poiList_;
float MultiDimFit::deltaNLL_ = 0;
unsigned int MultiDimFit::g_points = 50;
unsigned int MultiDimFit::firstPoint_ = 0;
unsigned int MultiDimFit::lastPoint_ = std::numeric_limits<unsigned int>::max();
bool MultiDimFit::floatOtherPOIs_ = false;
unsigned int MultiDimFit::nOtherFloatingPoi_ = 0;
bool MultiDimFit::fastScan_ = false;
bool MultiDimFit::loadedSnapshot_ = false;
bool MultiDimFit::savingSnapshot_ = false;
bool MultiDimFit::startFromPreFit_ = false;
bool MultiDimFit::alignEdges_ = false;
bool MultiDimFit::hasMaxDeltaNLLForProf_ = false;
bool MultiDimFit::squareDistPoiStep_ = false;
bool MultiDimFit::skipInitialFit_ = false;
bool MultiDimFit::saveFitResult_ = false;
float MultiDimFit::maxDeltaNLLForProf_ = 200;
float MultiDimFit::autoRange_ = -1.0;
std::string MultiDimFit::fixedPointPOIs_ = "";
float MultiDimFit::centeredRange_ = -1.0;
bool MultiDimFit::robustHesse_ = false;
std::string MultiDimFit::robustHesseLoad_ = "";
std::string MultiDimFit::robustHesseSave_ = "";

std::string MultiDimFit::saveSpecifiedFuncs_;
std::string MultiDimFit::saveSpecifiedIndex_;
std::string MultiDimFit::saveSpecifiedNuis_;
std::vector<std::string> MultiDimFit::specifiedFuncNames_;
std::vector<RooAbsReal *> MultiDimFit::specifiedFunc_;
std::vector<float> MultiDimFit::specifiedFuncVals_;
RooArgList MultiDimFit::specifiedFuncList_;
std::vector<std::string> MultiDimFit::specifiedCatNames_;
std::vector<RooCategory *> MultiDimFit::specifiedCat_;
std::vector<int> MultiDimFit::specifiedCatVals_;
RooArgList MultiDimFit::specifiedCatList_;
std::vector<std::string> MultiDimFit::specifiedNuis_;
std::vector<RooRealVar *> MultiDimFit::specifiedVars_;
std::vector<float> MultiDimFit::specifiedVals_;
RooArgList MultiDimFit::specifiedList_;
bool MultiDimFit::saveInactivePOI_ = false;

MultiDimFit::MultiDimFit() : FitterAlgoBase("MultiDimFit specific options") {
}

void MultiDimFit::applyOptions() {
  applyOptionsBase();
  g_algo = Algo::Grid;
}

bool MultiDimFit::runSpecific(RooWorkspace *w,
                              RooStats::ModelConfig *mc_s,
                              RooStats::ModelConfig *mc_b,
                              RooAbsData &data,
                              double &limit,
                              double &limitErr,
                              const double *hint) {
  // one-time initialization of POI variables, TTree branches, ...
  Combine::toggleGlobalFillTree(true);

  static int isInit = false;
  if (!isInit) {
    initOnce(w, mc_s);
    isInit = true;
  }

  // Get PDF
  RooAbsPdf &pdf = *mc_s->GetPdf();

  // Process POI not in list
  nOtherFloatingPoi_ = 0;
  int nConstPoi = 0;
  RooLinkedListIter iterP = mc_s->GetParametersOfInterest()->iterator();
  std::string setConstPOI;
  for (RooAbsArg *a = (RooAbsArg *)iterP.Next(); a != 0; a = (RooAbsArg *)iterP.Next()) {
    if (poiList_.contains(*a))
      continue;
    RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
    if (rrv == 0) {
      std::cerr << "combine/MultiDimFit: Parameter of interest " << a->GetName()
                << " which is not a RooRealVar will be ignored" << std::endl;
      continue;
    }
    rrv->setConstant(!floatOtherPOIs_);
    if (!floatOtherPOIs_) {
      setConstPOI += std::string(rrv->GetName()) + ", ";
      nConstPoi++;
    }
    if (floatOtherPOIs_)
      nOtherFloatingPoi_++;
  }
  if (nConstPoi > 0)
    std::cout << "Following POIs have been set constant (use --floatOtherPOIs to let them float): " << setConstPOI
              << std::endl;

  // start with a best fit
  const RooCmdArg &constrainCmdArg =
      g_withSystematics ? RooFit::Constrain(*mc_s->GetNuisanceParameters()) : RooCmdArg();
  std::unique_ptr<RooFitResult> res;
  if (g_verbose <= 3)
    RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CountErrors);
  bool doHesse = (g_algo == Algo::Singles || g_algo == Algo::Impact) || (saveFitResult_);
  if (!skipInitialFit_) {
    res.reset(doFit(pdf,
                    data,
                    (doHesse ? poiList_ : RooArgList()),
                    constrainCmdArg,
                    (saveFitResult_ && !robustHesse_),
                    1,
                    true,
                    false));
    if (!res.get()) {
      std::cout << "\n " << std::endl;
      std::cout << "\n ---------------------------" << std::endl;
      std::cout << "\n WARNING: MultiDimFit failed" << std::endl;
      std::cout << "\n ---------------------------" << std::endl;
      std::cout << "\n " << std::endl;
    }
    if (g_algo == Algo::Impact && res.get()) {
      // Set the floating parameters back to the best-fit value
      // before we write an entry into the output TTree
      w->allVars().assignValueOnly(res.get()->floatParsFinal());
    }
  } else {
    std::cout << "combine/MultiDimFit -- Skipping initial global fit" << std::endl;
    // must still create the NLL
    nll.reset(pdf.createNLL(data, constrainCmdArg, RooFit::Extended(pdf.canBeExtended()), RooFit::Offset(true)));
  }

  //if(w->var("r")) {w->var("r")->Print();}
  if (loadedSnapshot_ || res.get() || keepFailures_) {
    for (int i = 0, n = poi_.size(); i < n; ++i) {
      if (res.get() && doHesse) {
        // (res.get())->Print("v");
        RooAbsArg *rfloat = (res.get())->floatParsFinal().find(poi_[i].c_str());
        if (!rfloat) {
          rfloat = (*res).constPars().find(poi_[i].c_str());
        }
        RooRealVar *rf = dynamic_cast<RooRealVar *>(rfloat);
        poiVals_[i] = rf->getVal();  //for Singles we store the RooFitResults values
      } else
        poiVals_[i] = poiVars_[i]->getVal();
    }
    //if (g_algo != Algo::None) {
    for (unsigned int j = 0; j < specifiedNuis_.size(); j++) {
      specifiedVals_[j] = specifiedVars_[j]->getVal();
    }
    for (unsigned int j = 0; j < specifiedFuncNames_.size(); j++) {
      specifiedFuncVals_[j] = specifiedFunc_[j]->getVal();
    }
    for (unsigned int j = 0; j < specifiedCatNames_.size(); j++) {
      specifiedCatVals_[j] = specifiedCat_[j]->getIndex();
    }
    Combine::commitPoint(/*expected=*/false,
                         /*quantile=*/-1.);  // Combine will not commit a point anymore at -1 so can do it here
                                             //}
  }

  if (robustHesse_) {
    RobustHesse robustHesse(*nll, g_verbose - 1);
    robustHesse.ProtectArgSet(*mc_s->GetParametersOfInterest());
    if (robustHesseSave_ != "") {
      robustHesse.SaveHessianToFile(robustHesseSave_);
    }
    if (robustHesseLoad_ != "") {
      robustHesse.LoadHessianFromFile(robustHesseLoad_);
    }
    robustHesse.hesse();
    if (saveFitResult_) {
      res.reset(robustHesse.GetRooFitResult(res.get()));
    }
    robustHesse.WriteOutputFile("robustHesse" + name_ + ".root");
  }

  //set snapshot for best fit
  if (savingSnapshot_)
    w->saveSnapshot("MultiDimFit", utils::returnAllVars(w));

  if (autoRange_ > 0) {
    std::cout << "Adjusting range of POIs to +/- " << autoRange_ << " standard deviations" << std::endl;
    for (int i = 0, n = poi_.size(); i < n; ++i) {
      double val = poiVars_[i]->getVal(), err = poiVars_[i]->getError(), min0 = poiVars_[i]->getMin(),
             max0 = poiVars_[i]->getMax();
      double min1 = std::max(min0, val - autoRange_ * err);
      double max1 = std::min(max0, val + autoRange_ * err);
      std::cout << poi_[i] << ": " << val << " +/- " << err << " [ " << min0 << " , " << max0 << " ] ==> [ " << min1
                << " , " << max1 << " ]" << std::endl;
      poiVars_[i]->setRange(min1, max1);
    }
  }
  if (centeredRange_ > 0) {
    std::cout << "Adjusting range of POIs to +/- " << centeredRange_ << std::endl;
    for (int i = 0, n = poi_.size(); i < n; ++i) {
      double val = poiVars_[i]->getVal(), min0 = poiVars_[i]->getMin(), max0 = poiVars_[i]->getMax();
      double min1 = std::max(min0, val - centeredRange_);
      double max1 = std::min(max0, val + centeredRange_);
      std::cout << poi_[i] << ": " << val << " [ " << min0 << " , " << max0 << " ] ==> [ " << min1 << " , " << max1
                << " ]" << std::endl;
      poiVars_[i]->setRange(min1, max1);
    }
  }

  switch (g_algo) {
    case Algo::None: {
      std::cout << "\n --- MultiDimFit ---" << std::endl;
      std::cout << "best fit parameter values: " << std::endl;
      int len = poi_[0].length();
      for (int i = 0, n = poi_.size(); i < n; ++i) {
        len = std::max<int>(len, poi_[i].length());
      }
      for (int i = 0, n = poi_.size(); i < n; ++i) {
        printf("   %*s :  %+8.3f\n", len, poi_[i].c_str(), poiVals_[i]);
      }
    }
      if (res.get() && saveFitResult_)
        saveResult(*res);
      break;
    case Algo::Singles:
      if (res.get()) {
        doSingles(*res);
        if (saveFitResult_) {
          saveResult(*res);
        }
      }
      break;
    case Algo::Cross:
      doBox(*nll, g_confidenceLevel, "box", true);
      break;
    case Algo::Grid:
      doGrid(w, *nll);
      break;
    case Algo::RandomPoints:
      doRandomPoints(w, *nll);
      break;
    case Algo::FixedPoint:
      doFixedPoint(w, *nll);
      break;
    case Algo::Contour2D:
      doContour2D(w, *nll);
      break;
    case Algo::Stitch2D:
      doStitch2D(w, *nll);
      break;
    case Algo::Impact:
      if (res.get())
        doImpact(*res, *nll);
      break;
  }

  Combine::toggleGlobalFillTree(false);
  return true;
}

void MultiDimFit::initOnce(RooWorkspace *w, RooStats::ModelConfig *mc_s) {
  // Tell combine not to Fill its tree, we'll do it here;

  RooArgSet mcPoi(*mc_s->GetParametersOfInterest());
  if (poi_.empty()) {
    RooLinkedListIter iterP = mc_s->GetParametersOfInterest()->iterator();
    for (RooAbsArg *a = (RooAbsArg *)iterP.Next(); a != 0; a = (RooAbsArg *)iterP.Next()) {
      poi_.push_back(a->GetName());
    }
  }
  for (std::vector<std::string>::const_iterator it = poi_.begin(), ed = poi_.end(); it != ed; ++it) {
    RooAbsArg *a = mcPoi.find(it->c_str());
    bool isPoi = true;
    if (a == 0) {
      a = w->arg(it->c_str());  // look for the parameter elsewhere, but remember to clear its optimizeBounds attribute
      isPoi = false;
    }
    if (a == 0)
      throw std::invalid_argument(std::string("Parameter of interest ") + *it + " not in model.");
    RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
    if (rrv == 0)
      throw std::invalid_argument(std::string("Parameter of interest ") + *it + " not a RooRealVar.");
    if (!isPoi) {
      if (rrv->getAttribute("optimizeBounds")) {
        rrv->setAttribute("optimizeBounds", false);
        if (rrv->hasRange("optimizeBoundRange"))
          rrv->setRange(rrv->getMin("optimizeBoundRange"), rrv->getMax("optimizeBoundRange"));
      }
    }
    poiVars_.push_back(rrv);
    poiVals_.push_back(rrv->getVal());
    poiList_.add(*rrv);
  }

  if (saveSpecifiedFuncs_ != "") {
    char tmp[10240];
    strlcpy(tmp, saveSpecifiedFuncs_.c_str(), 10240);
    char *token = strtok(tmp, ",");
    while (token) {
      RooAbsArg *a = w->arg(token);
      if (a == 0)
        throw std::invalid_argument(std::string("function ") + token + " not in model.");
      RooAbsReal *rrv = dynamic_cast<RooAbsReal *>(a);
      if (rrv == 0)
        throw std::invalid_argument(std::string("function ") + token + " not a RooAbsReal.");
      specifiedFuncNames_.push_back(token);
      specifiedFunc_.push_back(rrv);
      specifiedFuncVals_.push_back(rrv->getVal());
      specifiedFuncList_.add(*rrv);
      token = strtok(0, ",");
    }
  }
  if (saveSpecifiedIndex_ != "") {
    char tmp[10240];
    strlcpy(tmp, saveSpecifiedIndex_.c_str(), 10240);
    char *token = strtok(tmp, ",");
    while (token) {
      RooCategory *rrv = w->cat(token);
      if (rrv == 0)
        throw std::invalid_argument(std::string("function ") + token + " not a RooCategory.");
      specifiedCatNames_.push_back(token);
      specifiedCat_.push_back(rrv);
      specifiedCatVals_.push_back(rrv->getIndex());
      specifiedCatList_.add(*rrv);
      token = strtok(0, ",");
    }
  }

  if (saveSpecifiedNuis_ != "" && g_withSystematics) {
    RooArgSet mcNuis(*mc_s->GetNuisanceParameters());
    if (saveSpecifiedNuis_ == "all") {
      specifiedNuis_.clear();
      RooLinkedListIter iterN = mc_s->GetNuisanceParameters()->iterator();
      for (RooAbsArg *a = (RooAbsArg *)iterN.Next(); a != 0; a = (RooAbsArg *)iterN.Next()) {
        specifiedNuis_.push_back(a->GetName());
      }
    } else {
      char tmp[10240];
      strlcpy(tmp, saveSpecifiedNuis_.c_str(), 10240);
      char *token = strtok(tmp, ",");
      while (token) {
        const RooArgSet *group = mc_s->GetWS()->set((std::string("group_") + token).data());
        if (group) {
          RooLinkedListIter iterN = group->iterator();
          for (RooAbsArg *a = (RooAbsArg *)iterN.Next(); a != 0; a = (RooAbsArg *)iterN.Next()) {
            specifiedNuis_.push_back(a->GetName());
          }
        } else {
          specifiedNuis_.push_back(token);
        }
        token = strtok(0, ",");
      }
    }
    for (std::vector<std::string>::const_iterator it = specifiedNuis_.begin(), ed = specifiedNuis_.end(); it != ed;
         ++it) {
      RooAbsArg *a = mcNuis.find(it->c_str());
      if (a == 0)
        throw std::invalid_argument(std::string("Nuisance Parameter ") + *it + " not in model.");
      if (poiList_.contains(*a))
        continue;
      RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
      if (rrv == 0)
        throw std::invalid_argument(std::string("Nuisance Parameter ") + *it + " not a RooRealVar.");
      specifiedVars_.push_back(rrv);
      specifiedVals_.push_back(rrv->getVal());
      specifiedList_.add(*rrv);
    }
  }
  if (saveInactivePOI_) {
    RooLinkedListIter iterP = mc_s->GetParametersOfInterest()->iterator();
    for (RooAbsArg *a = (RooAbsArg *)iterP.Next(); a != 0; a = (RooAbsArg *)iterP.Next()) {
      if (poiList_.contains(*a))
        continue;
      if (specifiedList_.contains(*a))
        continue;
      RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
      specifiedNuis_.push_back(a->GetName());
      specifiedVars_.push_back(rrv);
      specifiedVals_.push_back(rrv->getVal());
      specifiedList_.add(*rrv);
    }
  }

  // then add the branches to the tree (at the end, so there are no resizes)
  for (int i = 0, n = poi_.size(); i < n; ++i) {
    Combine::addBranch(poi_[i].c_str(), &poiVals_[i], (poi_[i] + "/F").c_str());
  }
  for (int i = 0, n = specifiedNuis_.size(); i < n; ++i) {
    Combine::addBranch(specifiedNuis_[i].c_str(), &specifiedVals_[i], (specifiedNuis_[i] + "/F").c_str());
  }
  for (int i = 0, n = specifiedFuncNames_.size(); i < n; ++i) {
    Combine::addBranch(specifiedFuncNames_[i].c_str(), &specifiedFuncVals_[i], (specifiedFuncNames_[i] + "/F").c_str());
  }
  for (int i = 0, n = specifiedCatNames_.size(); i < n; ++i) {
    Combine::addBranch(specifiedCatNames_[i].c_str(), &specifiedCatVals_[i], (specifiedCatNames_[i] + "/I").c_str());
  }
  Combine::addBranch("deltaNLL", &deltaNLL_, "deltaNLL/F");
}

void MultiDimFit::doSingles(RooFitResult &res) {
  std::cout << "\n --- MultiDimFit ---" << std::endl;
  std::cout << "best fit parameter values and profile-likelihood uncertainties: " << std::endl;
  int len = poi_[0].length();
  for (int i = 0, n = poi_.size(); i < n; ++i) {
    len = std::max<int>(len, poi_[i].length());
  }
  for (int i = 0, n = poi_.size(); i < n; ++i) {
    RooAbsArg *rfloat = res.floatParsFinal().find(poi_[i].c_str());
    if (!rfloat) {
      rfloat = res.constPars().find(poi_[i].c_str());
    }
    RooRealVar *rf = dynamic_cast<RooRealVar *>(rfloat);
    double bestFitVal = rf->getVal();

    double hiErr = +(rf->hasRange("err68") ? rf->getMax("err68") - bestFitVal : rf->getAsymErrorHi());
    double loErr = -(rf->hasRange("err68") ? rf->getMin("err68") - bestFitVal : rf->getAsymErrorLo());
    double maxError = std::max<double>(std::max<double>(hiErr, loErr), rf->getError());
    if (fabs(hiErr) < 0.001 * maxError) {
      std::cout << " Warning - No valid high-error found, will report difference to maximum of range for : "
                << rf->GetName() << std::endl;
      hiErr = -bestFitVal + rf->getMax();
    }
    if (fabs(loErr) < 0.001 * maxError) {
      std::cout << " Warning - No valid low-error found, will report difference to minimum of range for : "
                << rf->GetName() << std::endl;
      loErr = +bestFitVal - rf->getMin();
    }

    poiVals_[i] = bestFitVal - loErr;
    Combine::commitPoint(true, /*quantile=*/-0.32);
    poiVals_[i] = bestFitVal + hiErr;
    Combine::commitPoint(true, /*quantile=*/0.32);

    double hiErr95 = +(do95_ && rf->hasRange("err95") ? rf->getMax("err95") - bestFitVal : 0);
    double loErr95 = -(do95_ && rf->hasRange("err95") ? rf->getMin("err95") - bestFitVal : 0);
    double maxError95 = std::max<double>(std::max<double>(hiErr95, loErr95), 2 * rf->getError());

    if (do95_ && rf->hasRange("err95")) {
      if (fabs(hiErr95) < 0.001 * maxError95) {
        std::cout << " Warning - No valid high-error (for 95%) found, will report difference to maximum of range for : "
                  << rf->GetName() << std::endl;
        hiErr95 = -bestFitVal + rf->getMax();
      }
      if (fabs(loErr95) < 0.001 * maxError95) {
        std::cout << " Warning - No valid low-error (for 95%) found, will report difference to minimum of range for : "
                  << rf->GetName() << std::endl;
        loErr95 = +bestFitVal - rf->getMin();
      }
      poiVals_[i] = bestFitVal - loErr95;
      Combine::commitPoint(true, /*quantile=*/-0.05);
      poiVals_[i] = bestFitVal + hiErr95;
      Combine::commitPoint(true, /*quantile=*/0.05);
      //poiVals_[i] = rf->getMax("err95"); Combine::commitPoint(true, /*quantile=*/-0.05);
      //poiVals_[i] = rf->getMin("err95"); Combine::commitPoint(true, /*quantile=*/0.05);
      poiVals_[i] = bestFitVal;
      printf("   %*s :  %+8.3f   %+6.3f/%+6.3f (68%%)    %+6.3f/%+6.3f (95%%) \n",
             len,
             poi_[i].c_str(),
             poiVals_[i],
             -loErr,
             hiErr,
             loErr95,
             -hiErr95);
    } else {
      poiVals_[i] = bestFitVal;
      printf("   %*s :  %+8.3f   %+6.3f/%+6.3f (68%%)\n", len, poi_[i].c_str(), poiVals_[i], -loErr, hiErr);
    }
  }
}

void MultiDimFit::doImpact(RooFitResult &res, RooAbsReal &nll) {
  std::cout << "\n --- MultiDimFit ---" << std::endl;
  std::cout << "Parameter impacts: " << std::endl;

  // Save the initial parameters here to reset between NPs
  std::unique_ptr<RooArgSet> params(nll.getParameters((const RooArgSet *)0));
  RooArgSet init_snap;
  params->snapshot(init_snap);

  // Save the best-fit values of the saved parameters
  // we want to measure the impacts on
  std::vector<float> specifiedVals = specifiedVals_;
  std::vector<float> impactLo = specifiedVals_;
  std::vector<float> impactHi = specifiedVals_;

  int len = 9;
  for (int i = 0, n = poi_.size(); i < n; ++i) {
    len = std::max<int>(len, poi_[i].length());
  }
  printf("  %-*s :   %-21s", len, "Parameter", "Best-fit");
  for (int i = 0, n = specifiedNuis_.size(); i < n; ++i) {
    printf("  %-13s", specifiedNuis_[i].c_str());
  }
  printf("\n");

  for (int i = 0, n = poi_.size(); i < n; ++i) {
    RooAbsArg *rfloat = res.floatParsFinal().find(poi_[i].c_str());
    if (!rfloat) {
      rfloat = res.constPars().find(poi_[i].c_str());
    }
    RooRealVar *rf = dynamic_cast<RooRealVar *>(rfloat);
    double bestFitVal = rf->getVal();

    double hiErr = +(rf->hasRange("err68") ? rf->getMax("err68") - bestFitVal : rf->getAsymErrorHi());
    double loErr = -(rf->hasRange("err68") ? rf->getMin("err68") - bestFitVal : rf->getAsymErrorLo());
    printf("  %-*s : %+8.3f  %+6.3f/%+6.3f", len, poi_[i].c_str(), bestFitVal, -loErr, hiErr);
    // Reset all parameters to initial state
    *params = init_snap;
    // Then set this NP constant
    poiVars_[i]->setConstant(true);
    CascadeMinimizer minim(nll, CascadeMinimizer::Constrained);
    //minim.setStrategy(minimizerStrategy_);
    // Another snapshot to reset between high and low fits
    RooArgSet snap;
    params->snapshot(snap);
    std::vector<double> doVals = {bestFitVal - loErr, bestFitVal + hiErr};
    for (unsigned x = 0; x < doVals.size(); ++x) {
      *params = snap;
      poiVals_[i] = doVals[x];
      poiVars_[i]->setVal(doVals[x]);
      bool ok = minim.minimize(g_verbose - 1);
      if (ok) {
        for (unsigned int j = 0; j < poiVars_.size(); j++) {
          poiVals_[j] = poiVars_[j]->getVal();
        }
        for (unsigned int j = 0; j < specifiedNuis_.size(); j++) {
          specifiedVals_[j] = specifiedVars_[j]->getVal();
        }
        for (unsigned int j = 0; j < specifiedFuncNames_.size(); j++) {
          specifiedFuncVals_[j] = specifiedFunc_[j]->getVal();
        }
        for (unsigned int j = 0; j < specifiedCatNames_.size(); j++) {
          specifiedCatVals_[j] = specifiedCat_[j]->getIndex();
        }
        Combine::commitPoint(true, /*quantile=*/0.32);
      }
      for (unsigned int j = 0; j < specifiedNuis_.size(); j++) {
        if (x == 0) {
          impactLo[j] = specifiedVars_[j]->getVal() - specifiedVals[j];
        } else if (x == 1) {
          impactHi[j] = specifiedVars_[j]->getVal() - specifiedVals[j];
        }
      }
    }
    for (unsigned j = 0; j < specifiedVals.size(); ++j) {
      printf("  %+6.3f/%+6.3f", impactLo[j], impactHi[j]);
    }
    printf("\n");
  }
}

void MultiDimFit::doGrid(RooWorkspace *w, RooAbsReal &nll) {
  unsigned int n = poi_.size();
  //if (poi_.size() > 2) throw std::logic_error("Don't know how to do a grid with more than 2 POIs.");
  double nll0 = nll.getVal();
  if (startFromPreFit_)
    w->loadSnapshot("clean");

  std::vector<double> p0(n), pmin(n), pmax(n);
  for (unsigned int i = 0; i < n; ++i) {
    p0[i] = poiVars_[i]->getVal();
    pmin[i] = poiVars_[i]->getMin();
    pmax[i] = poiVars_[i]->getMax();
    poiVars_[i]->setConstant(true);
    std::cout << " POI: " << poiVars_[i]->GetName() << "= " << p0[i] << " -> [" << pmin[i] << "," << pmax[i] << "]"
              << std::endl;
  }

  CascadeMinimizer minim(nll, CascadeMinimizer::Constrained);
  if (!autoBoundsPOIs_.empty())
    minim.setAutoBounds(&autoBoundsPOISet_);
  if (!autoMaxPOIs_.empty())
    minim.setAutoMax(&autoMaxPOISet_);
  //minim.setStrategy(minimizerStrategy_);
  std::unique_ptr<RooArgSet> params(nll.getParameters((const RooArgSet *)0));
  RooArgSet snap;
  params->snapshot(snap);
  //snap.Print("V");
  if (n == 1) {
    double xspacing = (pmax[0] - pmin[0]) / g_points;
    double xspacingOffset = 0.5;
    if (alignEdges_) {
      xspacing = (pmax[0] - pmin[0]) / (g_points - 1);
      if (g_points == 1)
        xspacing = 0;
      xspacingOffset = 0.0;
    }
    // can do a more intellegent spacing of g_points
    double xbestpoint = (p0[0] - pmin[0]) / xspacing;
    if (lastPoint_ == std::numeric_limits<unsigned int>::max()) {
      lastPoint_ = g_points - 1;
    }
    for (unsigned int i = 0; i < g_points; ++i) {
      if (i < firstPoint_)
        continue;
      if (i > lastPoint_)
        break;
      double x = pmin[0] + (i + xspacingOffset) * xspacing;
      // If we're aligning with the edges and this is the last point,
      // set x to pmax[0] exactly
      if (alignEdges_ && i == (g_points - 1)) {
        x = pmax[0];
      }
      if (xbestpoint > lastPoint_) {
        int ireverse = lastPoint_ - i + firstPoint_;
        x = pmin[0] + (ireverse + xspacingOffset) * xspacing;
      }

      if (squareDistPoiStep_) {
        // distance between steps goes as ~square of distance from middle or range (could this be changed to from best fit value?)
        double phalf = (pmax[0] - pmin[0]) / 2;
        if (x < (pmin[0] + phalf)) {
          x = pmin[0] + TMath::Sqrt((x - pmin[0]) / phalf) * phalf;
        } else {
          x = pmax[0] - TMath::Sqrt((pmax[0] - x) / phalf) * phalf;
        }
      }

      //if (g_verbose > 1) std::cout << "Point " << i << "/" << g_points << " " << poiVars_[0]->GetName() << " = " << x << std::endl;
      std::cout << "Point " << i << "/" << g_points << " " << poiVars_[0]->GetName() << " = " << x << std::endl;
      *params = snap;
      poiVals_[0] = x;
      poiVars_[0]->setVal(x);
      // now we minimize
      nll.clearEvalErrorLog();
      deltaNLL_ = nll.getVal() - nll0;
      if (nll.numEvalErrors() > 0) {
        deltaNLL_ = 9990;
        for (unsigned int j = 0; j < specifiedNuis_.size(); j++) {
          specifiedVals_[j] = specifiedVars_[j]->getVal();
        }
        for (unsigned int j = 0; j < specifiedFuncNames_.size(); j++) {
          specifiedFuncVals_[j] = specifiedFunc_[j]->getVal();
        }
        Combine::commitPoint(true, /*quantile=*/0);
        continue;
      }
      bool ok = fastScan_ || (hasMaxDeltaNLLForProf_ && (nll.getVal() - nll0) > maxDeltaNLLForProf_)
                    ? true
                    : minim.minimize(g_verbose - 1);
      if (ok) {
        deltaNLL_ = nll.getVal() - nll0;
        double qN = 2 * (deltaNLL_);
        double prob = ROOT::Math::chisquared_cdf_c(qN, n + nOtherFloatingPoi_);
        for (unsigned int j = 0; j < specifiedNuis_.size(); j++) {
          specifiedVals_[j] = specifiedVars_[j]->getVal();
        }
        for (unsigned int j = 0; j < specifiedFuncNames_.size(); j++) {
          specifiedFuncVals_[j] = specifiedFunc_[j]->getVal();
        }
        for (unsigned int j = 0; j < specifiedCatNames_.size(); j++) {
          specifiedCatVals_[j] = specifiedCat_[j]->getIndex();
        }
        Combine::commitPoint(true, /*quantile=*/prob);
      }
    }
  } else if (n == 2) {
    unsigned int sqrn = ceil(sqrt(double(g_points)));
    unsigned int ipoint = 0, nprint = ceil(0.005 * sqrn * sqrn);
    RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CountErrors);
    CloseCoutSentry sentry(g_verbose < 2);
    double deltaX = (pmax[0] - pmin[0]) / sqrn;
    double deltaY = (pmax[1] - pmin[1]) / sqrn;
    double spacingOffset = 0.5;
    if (alignEdges_) {
      deltaX = (pmax[0] - pmin[0]) / (sqrn - 1);
      deltaY = (pmax[1] - pmin[1]) / (sqrn - 1);
      spacingOffset = 0.0;
      if (sqrn == 1) {
        deltaX = 0;
        deltaY = 0;
      }
    }
    for (unsigned int i = 0; i < sqrn; ++i) {
      for (unsigned int j = 0; j < sqrn; ++j, ++ipoint) {
        if (ipoint < firstPoint_)
          continue;
        if (ipoint > lastPoint_)
          break;
        *params = snap;
        double x = pmin[0] + (i + spacingOffset) * deltaX;
        double y = pmin[1] + (j + spacingOffset) * deltaY;
        if (g_verbose && (ipoint % nprint == 0)) {
          fprintf(sentry.trueStdOut(),
                  "Point %d/%d, (i,j) = (%d,%d), %s = %f, %s = %f\n",
                  ipoint,
                  sqrn * sqrn,
                  i,
                  j,
                  poiVars_[0]->GetName(),
                  x,
                  poiVars_[1]->GetName(),
                  y);
        }
        poiVals_[0] = x;
        poiVals_[1] = y;
        poiVars_[0]->setVal(x);
        poiVars_[1]->setVal(y);
        nll.clearEvalErrorLog();
        nll.getVal();
        if (nll.numEvalErrors() > 0) {
          for (unsigned int j = 0; j < specifiedNuis_.size(); j++) {
            specifiedVals_[j] = specifiedVars_[j]->getVal();
          }
          for (unsigned int j = 0; j < specifiedFuncNames_.size(); j++) {
            specifiedFuncVals_[j] = specifiedFunc_[j]->getVal();
          }
          for (unsigned int j = 0; j < specifiedCatNames_.size(); j++) {
            specifiedCatVals_[j] = specifiedCat_[j]->getIndex();
          }
          deltaNLL_ = 9999;
          Combine::commitPoint(true, /*quantile=*/0);
          if (gridType_ == G3x3) {
            for (int i2 = -1; i2 <= +1; ++i2) {
              for (int j2 = -1; j2 <= +1; ++j2) {
                if (i2 == 0 && j2 == 0)
                  continue;
                poiVals_[0] = x + 0.33333333 * i2 * deltaX;
                poiVals_[1] = y + 0.33333333 * j2 * deltaY;
                for (unsigned int j = 0; j < specifiedNuis_.size(); j++) {
                  specifiedVals_[j] = specifiedVars_[j]->getVal();
                }
                for (unsigned int j = 0; j < specifiedFuncNames_.size(); j++) {
                  specifiedFuncVals_[j] = specifiedFunc_[j]->getVal();
                }
                for (unsigned int j = 0; j < specifiedCatNames_.size(); j++) {
                  specifiedCatVals_[j] = specifiedCat_[j]->getIndex();
                }
                deltaNLL_ = 9999;
                Combine::commitPoint(true, /*quantile=*/0);
              }
            }
          }
          continue;
        }
        // now we minimize
        bool skipme = hasMaxDeltaNLLForProf_ && (nll.getVal() - nll0) > maxDeltaNLLForProf_;
        bool ok = fastScan_ || skipme ? true : minim.minimize(g_verbose - 1);
        if (ok) {
          deltaNLL_ = nll.getVal() - nll0;
          double qN = 2 * (deltaNLL_);
          double prob = ROOT::Math::chisquared_cdf_c(qN, n + nOtherFloatingPoi_);
          for (unsigned int j = 0; j < specifiedNuis_.size(); j++) {
            specifiedVals_[j] = specifiedVars_[j]->getVal();
          }
          for (unsigned int j = 0; j < specifiedFuncNames_.size(); j++) {
            specifiedFuncVals_[j] = specifiedFunc_[j]->getVal();
          }
          for (unsigned int j = 0; j < specifiedCatNames_.size(); j++) {
            specifiedCatVals_[j] = specifiedCat_[j]->getIndex();
          }
          Combine::commitPoint(true, /*quantile=*/prob);
        }
        if (gridType_ == G3x3) {
          bool forceProfile = !fastScan_ && std::min(fabs(deltaNLL_ - 1.15), fabs(deltaNLL_ - 2.995)) < 0.5;
          utils::CheapValueSnapshot center(*params);
          double x0 = x, y0 = y;
          for (int i2 = -1; i2 <= +1; ++i2) {
            for (int j2 = -1; j2 <= +1; ++j2) {
              if (i2 == 0 && j2 == 0)
                continue;
              center.writeTo(*params);
              x = x0 + 0.33333333 * i2 * deltaX;
              y = y0 + 0.33333333 * j2 * deltaY;
              poiVals_[0] = x;
              poiVars_[0]->setVal(x);
              poiVals_[1] = y;
              poiVars_[1]->setVal(y);
              nll.clearEvalErrorLog();
              nll.getVal();
              if (nll.numEvalErrors() > 0) {
                for (unsigned int j = 0; j < specifiedNuis_.size(); j++) {
                  specifiedVals_[j] = specifiedVars_[j]->getVal();
                }
                for (unsigned int j = 0; j < specifiedFuncNames_.size(); j++) {
                  specifiedFuncVals_[j] = specifiedFunc_[j]->getVal();
                }
                for (unsigned int j = 0; j < specifiedCatNames_.size(); j++) {
                  specifiedCatVals_[j] = specifiedCat_[j]->getIndex();
                }
                deltaNLL_ = 9999;
                Combine::commitPoint(true, /*quantile=*/0);
                continue;
              }
              deltaNLL_ = nll.getVal() - nll0;
              if (forceProfile || (!fastScan_ && std::min(fabs(deltaNLL_ - 1.15), fabs(deltaNLL_ - 2.995)) < 0.5)) {
                minim.minimize(g_verbose - 1);
                deltaNLL_ = nll.getVal() - nll0;
              }
              double qN = 2 * (deltaNLL_);
              double prob = ROOT::Math::chisquared_cdf_c(qN, n + nOtherFloatingPoi_);
              for (unsigned int j = 0; j < specifiedNuis_.size(); j++) {
                specifiedVals_[j] = specifiedVars_[j]->getVal();
              }
              for (unsigned int j = 0; j < specifiedFuncNames_.size(); j++) {
                specifiedFuncVals_[j] = specifiedFunc_[j]->getVal();
              }
              for (unsigned int j = 0; j < specifiedCatNames_.size(); j++) {
                specifiedCatVals_[j] = specifiedCat_[j]->getIndex();
              }
              Combine::commitPoint(true, /*quantile=*/prob);
            }
          }
        }
      }
    }

  } else {  // Use utils routine if n > 2

    unsigned int rootn = ceil(TMath::Power(double(g_points), double(1. / n)));
    unsigned int ipoint = 0, nprint = ceil(0.005 * TMath::Power((double)rootn, (double)n));

    RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CountErrors);
    CloseCoutSentry sentry(g_verbose < 2);

    // Create permutations
    std::vector<int> axis_g_points;

    for (unsigned int poi_i = 0; poi_i < n; poi_i++) {
      axis_g_points.push_back((int)rootn);
    }

    std::vector<std::vector<int> > permutations = utils::generateCombinations(axis_g_points);
    // Step through g_points
    std::vector<std::vector<int> >::iterator perm_it = permutations.begin();
    int npermutations = permutations.size();
    for (; perm_it != permutations.end(); perm_it++) {
      if (ipoint < firstPoint_) {
        ipoint++;
        continue;
      }
      if (ipoint > lastPoint_)
        break;
      *params = snap;

      if (g_verbose && (ipoint % nprint == 0)) {
        fprintf(sentry.trueStdOut(), "Point %d/%d, ", ipoint, npermutations);
      }
      for (unsigned int poi_i = 0; poi_i < n; poi_i++) {
        int ip = (*perm_it)[poi_i];
        double deltaXi = (pmax[poi_i] - pmin[poi_i]) / rootn;
        double spacingOffset = 0.5;
        if (alignEdges_) {
          deltaXi = (pmax[poi_i] - pmin[poi_i]) / (rootn - 1);
          if (rootn == 1) {
            deltaXi = 0.;
          }
          spacingOffset = 0.0;
        }
        double xi = pmin[poi_i] + deltaXi * (ip + spacingOffset);
        poiVals_[poi_i] = xi;
        poiVars_[poi_i]->setVal(xi);
        if (g_verbose && (ipoint % nprint == 0)) {
          fprintf(sentry.trueStdOut(), " %s = %f ", poiVars_[poi_i]->GetName(), xi);
        }
      }
      if (g_verbose && (ipoint % nprint == 0))
        fprintf(sentry.trueStdOut(), "\n");

      nll.clearEvalErrorLog();
      nll.getVal();
      if (nll.numEvalErrors() > 0) {
        for (unsigned int j = 0; j < specifiedNuis_.size(); j++) {
          specifiedVals_[j] = specifiedVars_[j]->getVal();
        }
        for (unsigned int j = 0; j < specifiedFuncNames_.size(); j++) {
          specifiedFuncVals_[j] = specifiedFunc_[j]->getVal();
        }
        for (unsigned int j = 0; j < specifiedCatNames_.size(); j++) {
          specifiedCatVals_[j] = specifiedCat_[j]->getIndex();
        }
        deltaNLL_ = 9999;
        Combine::commitPoint(true, /*quantile=*/0);
        ipoint++;
        continue;
      }
      // now we minimize
      bool skipme = hasMaxDeltaNLLForProf_ && (nll.getVal() - nll0) > maxDeltaNLLForProf_;
      bool ok = fastScan_ || skipme ? true : minim.minimize(g_verbose - 1);
      if (ok) {
        deltaNLL_ = nll.getVal() - nll0;
        double qN = 2 * (deltaNLL_);
        double prob = ROOT::Math::chisquared_cdf_c(qN, n + nOtherFloatingPoi_);
        for (unsigned int j = 0; j < specifiedNuis_.size(); j++) {
          specifiedVals_[j] = specifiedVars_[j]->getVal();
        }
        for (unsigned int j = 0; j < specifiedFuncNames_.size(); j++) {
          specifiedFuncVals_[j] = specifiedFunc_[j]->getVal();
        }
        for (unsigned int j = 0; j < specifiedCatNames_.size(); j++) {
          specifiedCatVals_[j] = specifiedCat_[j]->getIndex();
        }
        Combine::commitPoint(true, /*quantile=*/prob);
      }
      ipoint++;
    }
  }
}

void MultiDimFit::doRandomPoints(RooWorkspace *w, RooAbsReal &nll) {
  double nll0 = nll.getVal();
  if (startFromPreFit_)
    w->loadSnapshot("clean");
  for (unsigned int i = 0, n = poi_.size(); i < n; ++i) {
    poiVars_[i]->setConstant(true);
  }

  CascadeMinimizer minim(nll, CascadeMinimizer::Constrained);
  if (!autoBoundsPOIs_.empty())
    minim.setAutoBounds(&autoBoundsPOISet_);
  if (!autoMaxPOIs_.empty())
    minim.setAutoMax(&autoMaxPOISet_);
  //minim.setStrategy(minimizerStrategy_);
  unsigned int n = poi_.size();
  for (unsigned int j = 0; j < g_points; ++j) {
    for (unsigned int i = 0; i < n; ++i) {
      poiVars_[i]->randomize();
      poiVals_[i] = poiVars_[i]->getVal();
    }
    // now we minimize
    {
      CloseCoutSentry sentry(g_verbose < 3);
      bool ok = minim.minimize(g_verbose - 1);
      if (ok) {
        double qN = 2 * (nll.getVal() - nll0);
        double prob = ROOT::Math::chisquared_cdf_c(qN, n + nOtherFloatingPoi_);
        for (unsigned int j = 0; j < specifiedNuis_.size(); j++) {
          specifiedVals_[j] = specifiedVars_[j]->getVal();
        }
        for (unsigned int j = 0; j < specifiedFuncNames_.size(); j++) {
          specifiedFuncVals_[j] = specifiedFunc_[j]->getVal();
        }
        for (unsigned int j = 0; j < specifiedCatNames_.size(); j++) {
          specifiedCatVals_[j] = specifiedCat_[j]->getIndex();
        }
        Combine::commitPoint(true, /*quantile=*/prob);
      }
    }
  }
}
void MultiDimFit::doFixedPoint(RooWorkspace *w, RooAbsReal &nll) {
  double nll0 = nll.getVal();
  if (startFromPreFit_)
    w->loadSnapshot("clean");
  for (unsigned int i = 0, n = poi_.size(); i < n; ++i) {
    poiVars_[i]->setConstant(true);
  }

  CascadeMinimizer minim(nll, CascadeMinimizer::Constrained);
  if (!autoBoundsPOIs_.empty())
    minim.setAutoBounds(&autoBoundsPOISet_);
  if (!autoMaxPOIs_.empty())
    minim.setAutoMax(&autoMaxPOISet_);
  //minim.setStrategy(minimizerStrategy_);
  unsigned int n = poi_.size();

  //for (unsigned int i = 0; i < n; ++i) {
  //        std::cout<<" Before setting fixed point "<<poiVars_[i]->GetName()<<"= "<<poiVals_[i]<<std::endl;
  //}
  if (fixedPointPOIs_ != "") {
    utils::setModelParameters(fixedPointPOIs_, w->allVars());
  } else if (g_setPhysicsModelParameterExpression != "") {
    std::cout << " --fixedPointPOIs option not used, so will use the argument of --setParameters instead" << std::endl;
    utils::setModelParameters(g_setPhysicsModelParameterExpression, w->allVars());
  }

  for (unsigned int i = 0; i < n; ++i) {
    poiVars_[i]->setConstant(true);
    poiVals_[i] = poiVars_[i]->getVal();
    std::cout << " Evaluating fixed point with " << poiVars_[i]->GetName() << "= " << poiVals_[i] << std::endl;
  }
  // now we minimize
  {
    CloseCoutSentry sentry(g_verbose < 3);
    bool ok = minim.minimize(g_verbose - 1);
    if (ok) {
      nll0Value_ = nll0;
      nllValue_ = nll.getVal();
      deltaNLL_ = nll.getVal() - nll0;
      double qN = 2 * (nll.getVal() - nll0);
      double prob = ROOT::Math::chisquared_cdf_c(qN, n + nOtherFloatingPoi_);
      for (unsigned int j = 0; j < specifiedNuis_.size(); j++) {
        specifiedVals_[j] = specifiedVars_[j]->getVal();
      }
      for (unsigned int j = 0; j < specifiedFuncNames_.size(); j++) {
        specifiedFuncVals_[j] = specifiedFunc_[j]->getVal();
      }
      for (unsigned int j = 0; j < specifiedCatNames_.size(); j++) {
        specifiedCatVals_[j] = specifiedCat_[j]->getIndex();
      }
      Combine::commitPoint(true, /*quantile=*/prob);
      //for (unsigned int i = 0; i < n; ++i) {
      //        std::cout<<" after the fit "<<poiVars_[i]->GetName()<<"= "<<poiVars_[i]->getVal()<<std::endl;
      //}
    }
  }
}

void MultiDimFit::doContour2D(RooWorkspace *, RooAbsReal &nll) {
  if (poi_.size() != 2)
    throw std::logic_error("Contour2D works only in 2 dimensions");
  RooRealVar *xv = poiVars_[0];
  double x0 = poiVals_[0];
  float &x = poiVals_[0];
  RooRealVar *yv = poiVars_[1];
  double y0 = poiVals_[1];
  float &y = poiVals_[1];

  double threshold =
      nll.getVal() + 0.5 * ROOT::Math::chisquared_quantile_c(1 - g_confidenceLevel, 2 + nOtherFloatingPoi_);
  if (g_verbose > 0)
    std::cout << "Best fit point is for " << xv->GetName() << ", " << yv->GetName() << " =  " << x0 << ", " << y0
              << std::endl;

  // make a box
  doBox(nll, g_confidenceLevel, "box");
  double xMin = xv->getMin("box"), xMax = xv->getMax("box");
  double yMin = yv->getMin("box"), yMax = yv->getMax("box");

  g_verbose--;  // reduce verbosity to avoid messages from findCrossing
  // ===== Get relative min/max of x for several fixed y values =====
  yv->setConstant(true);
  for (unsigned int j = 0; j <= g_points; ++j) {
    if (j < firstPoint_)
      continue;
    if (j > lastPoint_)
      break;
    // take g_points uniformly spaced in polar angle in the case of a perfect circle
    double yc = 0.5 * (yMax + yMin), yr = 0.5 * (yMax - yMin);
    yv->setVal(yc + yr * std::cos(j * M_PI / double(g_points)));
    // ===== Get the best fit x (could also do without profiling??) =====
    xv->setConstant(false);
    xv->setVal(x0);
    CascadeMinimizer minimXI(nll, CascadeMinimizer::Unconstrained, xv);
    if (!autoBoundsPOIs_.empty())
      minimXI.setAutoBounds(&autoBoundsPOISet_);
    if (!autoMaxPOIs_.empty())
      minimXI.setAutoMax(&autoMaxPOISet_);
    //minimXI.setStrategy(minimizerStrategy_);
    {
      CloseCoutSentry sentry(g_verbose < 3);
      minimXI.minimize(g_verbose - 1);
    }
    double xc = xv->getVal();
    xv->setConstant(true);
    if (g_verbose > -1)
      std::cout << "Best fit " << xv->GetName() << " for  " << yv->GetName() << " = " << yv->getVal() << " is at " << xc
                << std::endl;
    // ===== Then get the range =====
    CascadeMinimizer minim(nll, CascadeMinimizer::Constrained);
    if (!autoBoundsPOIs_.empty())
      minim.setAutoBounds(&autoBoundsPOISet_);
    if (!autoMaxPOIs_.empty())
      minim.setAutoMax(&autoMaxPOISet_);
    double xup = findCrossing(minim, nll, *xv, threshold, xc, xMax);
    if (!std::isnan(xup)) {
      x = xup;
      y = yv->getVal();
      Combine::commitPoint(true, /*quantile=*/1 - g_confidenceLevel);
      if (g_verbose > -1)
        std::cout << "Minimum of " << xv->GetName() << " at " << g_confidenceLevel << " CL for " << yv->GetName()
                  << " = " << y << " is " << x << std::endl;
    }

    double xdn = findCrossing(minim, nll, *xv, threshold, xc, xMin);
    if (!std::isnan(xdn)) {
      x = xdn;
      y = yv->getVal();
      Combine::commitPoint(true, /*quantile=*/1 - g_confidenceLevel);
      if (g_verbose > -1)
        std::cout << "Maximum of " << xv->GetName() << " at " << g_confidenceLevel << " CL for " << yv->GetName()
                  << " = " << y << " is " << x << std::endl;
    }
  }

  g_verbose++;  // restore verbosity
}

void MultiDimFit::doStitch2D(RooWorkspace *, RooAbsReal &nll) {
  if (poi_.size() != 2)
    throw std::logic_error("Contour2D works only in 2 dimensions");
}

void MultiDimFit::doBox(RooAbsReal &nll, double g_confidenceLevel, const char *name, bool commitPoints) {
  unsigned int n = poi_.size();
  double nll0 = nll.getVal(),
         threshold = nll0 + 0.5 * ROOT::Math::chisquared_quantile_c(1 - g_confidenceLevel, n + nOtherFloatingPoi_);

  std::vector<double> p0(n);
  for (unsigned int i = 0; i < n; ++i) {
    p0[i] = poiVars_[i]->getVal();
    poiVars_[i]->setConstant(false);
  }

  g_verbose--;  // reduce verbosity due to findCrossing
  for (unsigned int i = 0; i < n; ++i) {
    RooRealVar *xv = poiVars_[i];
    xv->setConstant(true);
    CascadeMinimizer minimX(nll, CascadeMinimizer::Constrained);
    if (!autoBoundsPOIs_.empty())
      minimX.setAutoBounds(&autoBoundsPOISet_);
    if (!autoMaxPOIs_.empty())
      minimX.setAutoMax(&autoMaxPOISet_);
    //minimX.setStrategy(minimizerStrategy_);

    for (unsigned int j = 0; j < n; ++j)
      poiVars_[j]->setVal(p0[j]);
    double xMin = findCrossing(minimX, nll, *xv, threshold, p0[i], xv->getMin());
    if (!std::isnan(xMin)) {
      if (g_verbose > -1)
        std::cout << "Minimum of " << xv->GetName() << " at " << g_confidenceLevel << " CL for all others floating is "
                  << xMin << std::endl;
      for (unsigned int j = 0; j < n; ++j)
        poiVals_[j] = poiVars_[j]->getVal();
      if (commitPoints)
        Combine::commitPoint(true, /*quantile=*/1 - g_confidenceLevel);
    } else {
      xMin = xv->getMin();
      for (unsigned int j = 0; j < n; ++j)
        poiVals_[j] = poiVars_[j]->getVal();
      double prob = ROOT::Math::chisquared_cdf_c(2 * (nll.getVal() - nll0), n + nOtherFloatingPoi_);
      if (commitPoints)
        Combine::commitPoint(true, /*quantile=*/prob);
      if (g_verbose > -1)
        std::cout << "Minimum of " << xv->GetName() << " at " << g_confidenceLevel << " CL for all others floating is "
                  << xMin << " (on the boundary, p-val " << prob << ")" << std::endl;
    }

    for (unsigned int j = 0; j < n; ++j)
      poiVars_[j]->setVal(p0[j]);
    double xMax = findCrossing(minimX, nll, *xv, threshold, p0[i], xv->getMax());
    if (!std::isnan(xMax)) {
      if (g_verbose > -1)
        std::cout << "Maximum of " << xv->GetName() << " at " << g_confidenceLevel << " CL for all others floating is "
                  << xMax << std::endl;
      for (unsigned int j = 0; j < n; ++j)
        poiVals_[j] = poiVars_[j]->getVal();
      if (commitPoints)
        Combine::commitPoint(true, /*quantile=*/1 - g_confidenceLevel);
    } else {
      xMax = xv->getMax();
      double prob = ROOT::Math::chisquared_cdf_c(2 * (nll.getVal() - nll0), n + nOtherFloatingPoi_);
      for (unsigned int j = 0; j < n; ++j)
        poiVals_[j] = poiVars_[j]->getVal();
      if (commitPoints)
        Combine::commitPoint(true, /*quantile=*/prob);
      if (g_verbose > -1)
        std::cout << "Maximum of " << xv->GetName() << " at " << g_confidenceLevel << " CL for all others floating is "
                  << xMax << " (on the boundary, p-val " << prob << ")" << std::endl;
    }

    xv->setRange(name, xMin, xMax);
    xv->setConstant(false);
  }
  g_verbose++;  // restore verbosity
}

void MultiDimFit::saveResult(RooFitResult &res) {
  if (g_verbose > 2)
    res.Print("V");
  fitOut.reset(TFile::Open(("multidimfit" + name_ + ".root").c_str(), "RECREATE"));
  fitOut->WriteTObject(&res, "fit_mdf");
  fitOut->cd();
  fitOut.release()->Close();
}
