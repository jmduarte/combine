#include "combine/Significance.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooRandom.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
//#include "RooMinimizerOpt.h"
#include "RooMinimizer.h"
#include "TGraph.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TCanvas.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/HypoTestResult.h"
#include "RooStats/RooStatsUtils.h"
#include "combine/Combine.h"
#include "combine/CloseCoutSentry.h"
#include "combine/utils.h"
#include "combine/ProfiledLikelihoodRatioTestStatExt.h"
#include "combine/CascadeMinimizer.h"

#include <Math/MinimizerOptions.h>

using namespace RooStats;
using namespace std;

//std::string Significance::minimizerAlgo_ = "Minuit2";
std::string Significance::minimizerAlgoForBF_ = "Minuit2,simplex";
//float       Significance::minimizerTolerance_ = 1e-2;
float Significance::minimizerToleranceForBF_ = 1e-4;
int Significance::tries_ = 1;
int Significance::maxTries_ = 1;
float Significance::maxRelDeviation_ = 0.05;
float Significance::maxOutlierFraction_ = 0.25;
int Significance::maxOutliers_ = 3;
int Significance::points_ = 20;
bool Significance::preFit_ = false;
bool Significance::useMinos_ = true;
bool Significance::bruteForce_ = false;
std::string Significance::bfAlgo_ = "scale";
bool Significance::uncapped_ = false;
float Significance::signalForSignificance_ = 0;
std::string Significance::plot_ = "";

Significance::Significance() : LimitAlgo("Significance specific options") {
}

void Significance::applyOptions() {}

Significance::MinimizerSentry::MinimizerSentry(const std::string &minimizerAlgo, double tolerance)
    : minimizerTypeBackup(ROOT::Math::MinimizerOptions::DefaultMinimizerType()),
      minimizerAlgoBackup(ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo()),
      minimizerTollBackup(ROOT::Math::MinimizerOptions::DefaultTolerance()) {
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(tolerance);
  if (minimizerAlgo.find(",") != std::string::npos) {
    size_t idx = minimizerAlgo.find(",");
    std::string type = minimizerAlgo.substr(0, idx), algo = minimizerAlgo.substr(idx + 1);
    if (g_verbose > 1)
      std::cout << "Set default minimizer to " << type << ", algorithm " << algo << ", tolerance " << tolerance
                << std::endl;
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer(type.c_str(), algo.c_str());
  } else {
    if (g_verbose > 1)
      std::cout << "Set default minimizer to " << minimizerAlgo << ", tolerance " << tolerance << std::endl;
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer(minimizerAlgo.c_str());
  }
}

Significance::MinimizerSentry::~MinimizerSentry() {
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(minimizerTollBackup);
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer(minimizerTypeBackup.c_str(),
                                                    minimizerAlgoBackup.empty() ? 0 : minimizerAlgoBackup.c_str());
}

bool Significance::run(RooWorkspace *w,
                       RooStats::ModelConfig *mc_s,
                       RooStats::ModelConfig *mc_b,
                       RooAbsData &data,
                       double &limit,
                       double &limitErr,
                       const double *hint) {
  //MinimizerSentry minimizerConfig(minimizerAlgo_, minimizerTolerance_);
  CloseCoutSentry sentry(g_verbose < 0);

  double minimizerTolerance_ = ROOT::Math::MinimizerOptions::DefaultTolerance();
  RooRealVar *r = dynamic_cast<RooRealVar *>(mc_s->GetParametersOfInterest()->first());
  bool success = false;
  std::vector<double> limits;
  double rMax = r->getMax();
  std::unique_ptr<RooAbsPdf> nuisancePdf = nullptr;
  for (int i = 0; i < maxTries_; ++i) {
    w->loadSnapshot("clean");
    if (i > 0) {  // randomize starting point
      r->setMax(rMax * (0.5 + RooRandom::uniform()));
      r->setVal((0.1 + 0.5 * RooRandom::uniform()) * r->getMax());
      if (g_withSystematics) {
        if (nuisancePdf.get() == 0)
          nuisancePdf.reset(utils::makeNuisancePdf(*mc_s));
        RooArgSet set(*mc_s->GetNuisanceParameters());
        RooDataSet *randoms = nuisancePdf->generate(set, 1);
        set = *randoms->get(0);
        if (g_verbose > 2) {
          std::cout << "Starting minimization from point " << std::endl;
          r->Print("V");
          set.Print("V");
        }
        delete randoms;
      }
    }
    if (preFit_) {
      CloseCoutSentry sentry(g_verbose < 2);
      RooFitResult *res = mc_s->GetPdf()->fitTo(data, RooFit::Save(1), RooFit::Minimizer("Minuit2"));
      if (res == 0 || res->covQual() != 3 || res->edm() > minimizerTolerance_) {
        if (g_verbose > 1)
          std::cout << "Fit failed (covQual " << (res ? res->covQual() : -1) << ", edm " << (res ? res->edm() : 0)
                    << ")" << std::endl;
        continue;
      }
      if (g_verbose > 1) {
        res->Print("V");
        std::cout << "Covariance quality: " << res->covQual() << ", Edm = " << res->edm() << std::endl;
      }
      delete res;
    }
    //bool thisTry = (g_doSignificance ?  runSignificance(w,mc_s,data,limit,limitErr) : runLimit(w,mc_s,data,limit,limitErr)); <- only run signficance here now
    bool thisTry = runSignificance(w, mc_s, data, limit, limitErr);
    if (!thisTry)
      continue;
    if (tries_ == 1) {
      success = true;
      break;
    }
    limits.push_back(limit);
    int nresults = limits.size();
    if (nresults < tries_)
      continue;
    std::sort(limits.begin(), limits.end());
    double median = (nresults % 2 ? limits[nresults / 2] : 0.5 * (limits[nresults / 2] + limits[nresults / 2 + 1]));
    int noutlier = 0;
    double spreadIn = 0, spreadOut = 0;
    for (int j = 0; j < nresults; ++j) {
      double diff = fabs(limits[j] - median) / median;
      if (diff < maxRelDeviation_) {
        spreadIn = max(spreadIn, diff);
      } else {
        noutlier++;
        spreadOut = max(spreadOut, diff);
      }
    }
    if (g_verbose > 0) {
      std::cout << "Numer of tries: " << i << "   Number of successes: " << nresults << ", Outliers: " << noutlier
                << " (frac = " << noutlier / double(nresults) << ")"
                << ", Spread of non-outliers: " << spreadIn << " / of outliers: " << spreadOut << std::endl;
    }
    if (noutlier <= maxOutlierFraction_ * nresults) {
      if (g_verbose > 0)
        std::cout << " \\--> success! " << std::endl;
      success = true;
      limit = median;
      break;
    } else if (noutlier > maxOutliers_) {
      if (g_verbose > 0)
        std::cout << " \\--> failure! " << std::endl;
      break;
    }
  }
  return success;
}

bool Significance::runLimit(
    RooWorkspace *w, RooStats::ModelConfig *mc_s, RooAbsData &data, double &limit, double &limitErr) {
  RooArgSet poi(*mc_s->GetParametersOfInterest());
  RooRealVar *r = dynamic_cast<RooRealVar *>(poi.first());
  double rMax = r->getMax();
  bool success = false;
  CloseCoutSentry coutSentry(g_verbose <=
                             1);  // close standard output and error, so that we don't flood them with minuit messages

  double minimizerTolerance_ = ROOT::Math::MinimizerOptions::DefaultTolerance();

  while (!success) {
    ProfileLikelihoodCalculator plcB(data, *mc_s, 1.0 - g_confidenceLevel);
    std::unique_ptr<LikelihoodInterval> plInterval;
    if (useMinos_ || bruteForce_) {
      // try first with Minos, unless brute force requested
      if (!bruteForce_) {
        limit = upperLimitWithMinos(
            *mc_s->GetPdf(), data, *r, mc_s->GetNuisanceParameters(), minimizerTolerance_, g_confidenceLevel);
      }
      // if brute force forced, or minos failed, go with the next one
      if (std::isnan(limit) || bruteForce_) {
        std::pair<double, double> le = upperLimitBruteForce(
            *mc_s->GetPdf(), data, *r, mc_s->GetNuisanceParameters(), 1e-3 * minimizerTolerance_, g_confidenceLevel);
        limit = le.first;
        limitErr = le.second;
      }
    } else {
      plInterval.reset(plcB.GetInterval());
      if (plInterval.get() == 0)
        break;
      limit = g_lowerLimit ? plInterval->LowerLimit(*r) : plInterval->UpperLimit(*r);
    }
    if (limit >= 0.75 * r->getMax()) {
      std::cout << "Limit " << r->GetName() << " < " << limit << "; " << r->GetName() << " max < " << r->getMax()
                << std::endl;
      if (r->getMax() / rMax > 20)
        break;
      r->setMax(r->getMax() * 2);
      continue;
    }
    if (limit == r->getMin()) {
      std::cerr << "ProfileLikelihoodCalculator failed (returned upper limit equal to the lower bound)" << std::endl;
      break;
    }
    success = true;
    if (!plot_.empty()) {
      TCanvas *c1 = new TCanvas("c1", "c1");
      LikelihoodIntervalPlot plot(&*plInterval);
      plot.Draw();
      c1->Print(plot_.c_str());
      delete c1;
    }
  }
  coutSentry.clear();
  if (g_verbose >= 0) {
    if (success) {
      std::cout << "\n -- Significance -- "
                << "\n";
      if (limitErr) {
        std::cout << "Limit: " << r->GetName() << (g_lowerLimit ? " > " : " < ") << limit << " +/- " << limitErr
                  << " @ " << g_confidenceLevel * 100 << "% CL" << std::endl;
      } else {
        std::cout << "Limit: " << r->GetName() << (g_lowerLimit ? " > " : " < ") << limit << " @ "
                  << g_confidenceLevel * 100 << "% CL" << std::endl;
      }
    }
  }
  return success;
}

bool Significance::runSignificance(
    RooWorkspace *w, RooStats::ModelConfig *mc_s, RooAbsData &data, double &limit, double &limitErr) {
  RooArgSet poi(*mc_s->GetParametersOfInterest());
  RooRealVar *r = dynamic_cast<RooRealVar *>(poi.first());

  if (!uncapped_)
    r->setMin(signalForSignificance_);
  ProfileLikelihoodCalculator plcS(data, *mc_s, 1.0 - g_confidenceLevel);
  RooArgSet nullParamValues;
  r->setVal(signalForSignificance_);
  nullParamValues.addClone(*r);
  plcS.SetNullParameters(nullParamValues);

  CloseCoutSentry coutSentry(g_verbose <=
                             1);  // close standard output and error, so that we don't flood them with minuit messages
  double minimizerTolerance_ = ROOT::Math::MinimizerOptions::DefaultTolerance();

  if (bruteForce_) {
    double q0 = -1;
    if (bfAlgo_ == "scale")
      q0 = significanceBruteForce(*mc_s->GetPdf(), data, *r, mc_s->GetNuisanceParameters(), 0.1 * minimizerTolerance_);
    else
      q0 = significanceFromScan(
          *mc_s->GetPdf(), data, *r, mc_s->GetNuisanceParameters(), 0.1 * minimizerTolerance_, points_);
    if (q0 == -1)
      return false;
    limit = (q0 > 0 ? sqrt(2 * q0) : (uncapped_ ? -sqrt(-2 * q0) : 0));
  } else if (useMinos_) {
    ProfiledLikelihoodTestStatOpt testStat(
        *mc_s->GetObservables(),
        *mc_s->GetPdf(),
        mc_s->GetNuisanceParameters(),
        nullParamValues,
        *mc_s->GetParametersOfInterest(),
        RooArgList(),
        RooArgList(),
        g_verbose - 1,
        uncapped_ ? ProfiledLikelihoodTestStatOpt::signFlipDef : ProfiledLikelihoodTestStatOpt::oneSidedDef);
    Double_t q0 = testStat.Evaluate(data, nullParamValues);
    limit = (q0 > 0 ? sqrt(2 * q0) : (uncapped_ ? -sqrt(-2 * q0) : 0));
  } else {
    std::unique_ptr<HypoTestResult> result(plcS.GetHypoTest());
    if (result.get() == 0)
      return false;
    limit = result->Significance();
    if (uncapped_ && r->getVal() < signalForSignificance_)
      limit = -limit;
  }
  coutSentry.clear();
  if (limit == 0 && signbit(limit)) {
    //..... This is not an error, it just means we have a deficit of events.....
    std::cerr << "The minimum of the likelihood is for r <= " << signalForSignificance_
              << ", so the significance is zero" << std::endl;
    limit = 0;
  }
  return true;
}

double Significance::upperLimitWithMinos(RooAbsPdf &pdf,
                                         RooAbsData &data,
                                         RooRealVar &poi,
                                         const RooArgSet *nuisances,
                                         double tolerance,
                                         double confidenceLevel) const {
  std::unique_ptr<RooAbsReal> nll(pdf.createNLL(data, RooFit::Constrain(*nuisances)));
  RooMinimizer minim(*nll);
  minim.setStrategy(0);
  minim.setPrintLevel(g_verbose - 1);
  minim.setErrorLevel(0.5 * TMath::ChisquareQuantile(confidenceLevel, 1));
  nllutils::robustMinimize(*nll, minim, g_verbose - 1);
  int minosStat = minim.minos(RooArgSet(poi));
  if (minosStat == -1)
    return std::numeric_limits<double>::quiet_NaN();
  std::unique_ptr<RooFitResult> res(minim.save());
  double muhat = poi.getVal(), limit = poi.getVal() + (g_lowerLimit ? poi.getAsymErrorLo() : poi.getAsymErrorHi());
  double nll0 = nll->getVal();
  poi.setVal(limit);
  double nll2 = nll->getVal();
  if (nll2 < nll0 + 0.75 * 0.5 * TMath::ChisquareQuantile(confidenceLevel, 1)) {
    std::cerr << "ERROR: unprofiled likelihood gives better result than profiled one. deltaNLL = " << (nll2 - nll0)
              << ". will try brute force." << std::endl;
    poi.setVal(muhat);
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (g_verbose > 1)
    res->Print("V");
  return limit;
}

std::pair<double, double> Significance::upperLimitBruteForce(RooAbsPdf &pdf,
                                                             RooAbsData &data,
                                                             RooRealVar &poi,
                                                             const RooArgSet *nuisances,
                                                             double tolerance,
                                                             double confidenceLevel) const {
  poi.setConstant(false);
  std::unique_ptr<RooAbsReal> nll(pdf.createNLL(data, RooFit::Constrain(*nuisances)));
  RooMinimizer minim0(*nll);
  minim0.setStrategy(0);
  minim0.setPrintLevel(-1);
  nllutils::robustMinimize(*nll, minim0, g_verbose - 2);
  poi.setConstant(true);
  RooMinimizer minim(*nll);
  minim.setPrintLevel(-1);
  if (!nllutils::robustMinimize(*nll, minim, g_verbose - 2)) {
    std::cerr << "Initial minimization failed. Aborting." << std::endl;
    return std::pair<double, double>(0, -1);
  }
  std::unique_ptr<RooFitResult> start(minim.save());
  double minnll = nll->getVal();
  double rval = poi.getVal() + (g_lowerLimit ? -3 : +3) * poi.getError(), rlow = poi.getVal(),
         rhigh = g_lowerLimit ? poi.getMin() : poi.getMax();
  if (rval >= rhigh || rval <= rlow)
    rval = 0.5 * (rlow + rhigh);
  double target = minnll + 0.5 * TMath::ChisquareQuantile(confidenceLevel, 1);
  //minim.setPrintLevel(g_verbose-2);
  MinimizerSentry minimizerConfig(minimizerAlgoForBF_, minimizerToleranceForBF_);
  bool fail = false;
  if (g_verbose) {
    printf("  %-6s  delta(NLL)\n", poi.GetName());
    printf("%8.5f  %8.5f\n", rval, 0.);
    fflush(stdout);
  }
  do {
    poi.setVal(rval);
    minim.setStrategy(0);
    bool success = nllutils::robustMinimize(*nll, minim, g_verbose - 2);
    if (success == false) {
      std::cerr << "Minimization failed at " << poi.getVal() << ". exiting the bisection loop" << std::endl;
      fail = true;
      break;
    }
    double nllthis = nll->getVal();
    if (g_verbose) {
      printf("%8.5f  %8.5f\n", rval, nllthis - minnll);
      fflush(stdout);
    }
    if (fabs(nllthis - target) < tolerance) {
      return std::pair<double, double>(rval, (rhigh - rlow) * 0.5);
    } else if (nllthis < target) {
      (g_lowerLimit ? rhigh : rlow) = rval;
      rval = 0.5 * (rval + rhigh);
    } else {
      (g_lowerLimit ? rlow : rhigh) = rval;
      rval = 0.5 * (rval + rlow);
    }
  } while (fabs(rhigh - rlow) > tolerance);
  if (fail) {
    // try do do it in small steps instead
    std::unique_ptr<RooArgSet> pars(nll->getParameters((const RooArgSet *)0));
    double dx = (g_lowerLimit ? -0.05 : +0.05) * poi.getError();
    *pars = start->floatParsFinal();
    rval = poi.getVal() + dx;
    do {
      poi.setVal(rval);
      minim.setStrategy(0);
      bool success = nllutils::robustMinimize(*nll, minim, g_verbose - 2);
      if (success == false) {
        std::cerr << "Minimization failed at " << poi.getVal() << ". exiting the stepping loop" << std::endl;
        return std::pair<double, double>(poi.getVal(), fabs(rhigh - rlow) * 0.5);
      }
      double nllthis = nll->getVal();
      if (g_verbose) {
        printf("%8.5f  %8.5f\n", rval, nllthis - minnll);
        fflush(stdout);
      }
      if (fabs(nllthis - target) < tolerance) {
        return std::pair<double, double>(rval, fabs(dx));
      } else if (nllthis < target) {
        rval += dx;
      } else {
        dx *= 0.5;
        rval -= dx;
      }
    } while (rval < poi.getMax() && rval > poi.getMin());
    return std::pair<double, double>(poi.getMax(), 0);
  } else {
    return std::pair<double, double>(poi.getVal(), fabs(rhigh - rlow) * 0.5);
  }
}

double Significance::significanceBruteForce(
    RooAbsPdf &pdf, RooAbsData &data, RooRealVar &poi, const RooArgSet *nuisances, double tolerance) const {
  poi.setConstant(false);
  //poi.setMin(0);
  poi.setVal(0.05 * poi.getMax());
  std::unique_ptr<RooAbsReal> nll(pdf.createNLL(data, RooFit::Constrain(*nuisances)));
  CascadeMinimizer minim0(*nll, CascadeMinimizer::Unconstrained, &poi);
  minim0.setStrategy(0);
  minim0.minimize(g_verbose - 2);
  if (poi.getVal() < 0 && !uncapped_) {
    printf("Minimum found at %s = %8.5f < 0: significance will be zero\n", poi.GetName(), poi.getVal());
    return 0;
  }
  poi.setConstant(true);
  CascadeMinimizer minim(*nll, CascadeMinimizer::Constrained);
  if (!minim.minimize(g_verbose - 2)) {
    std::cerr << "Initial minimization failed. Aborting." << std::endl;
    return -1;
  } else if (g_verbose > 0) {
    printf("Minimum found at %s = %8.5f\n", poi.GetName(), poi.getVal());
  }
  MinimizerSentry minimizerConfig(minimizerAlgoForBF_, minimizerToleranceForBF_);
  std::unique_ptr<RooFitResult> start(minim.save());
  double minnll = nll->getVal(), thisnll = minnll, lastnll = thisnll;
  double rbest = poi.getVal(), rval = rbest;
  TGraph *points = 0;
  if (g_verbose) {
    printf("  %-6s  delta(NLL)\n", poi.GetName());
    printf("%8.5f  %8.5f\n", rval, 0.);
    fflush(stdout);
    points = new TGraph(1);
    points->SetName(Form("nll_scan_%g", g_mass));
    points->SetPoint(0, rval, 0);
  }
  while (std::abs(rval) >= tolerance * std::abs(rval > 0 ? poi.getMax() : poi.getMin())) {
    rval *= 0.8;
    poi.setVal(rval);
    minim.setStrategy(0);
    bool success = minim.improve(g_verbose - 2, /*cascade=*/false);
    lastnll = thisnll;
    thisnll = nll->getVal();
    if (success == false) {
      std::cerr << "Minimization failed at " << poi.getVal() << ". exiting the loop" << std::endl;
      return -1;
    }
    if (g_verbose) {
      printf("%8.5f  %8.5f\n", rval, thisnll - minnll);
      fflush(stdout);
      points->Set(points->GetN() + 1);
      points->SetPoint(points->GetN() - 1, rval, thisnll - minnll);
    }
    if (fabs(lastnll - thisnll) < 7 * minimizerToleranceForBF_) {
      std::cout << "This is enough." << std::endl;
      if (thisnll < lastnll) {
        std::cout << "Linear extrapolation from " << thisnll << " to " << (thisnll - (lastnll - thisnll) * rval)
                  << std::endl;
        thisnll -= (lastnll - thisnll) * rval;
      }
      break;
    }
#if 0
        if (lastnll > thisnll) {
            std::cout << "Secondary minimum found, back to original minimization" << std::endl;
            {
                MinimizerSentry minimizerConfig(minimizerAlgo_, minimizerTolerance_);
                poi.setConstant(false);
                poi.setVal(rbest);
                minim.improve(g_verbose - 1);
                printf("New global minimum at %8.5f (was %8.5f), the shift in nll is %8.5f\n", poi.getVal(), rbest, nll->getVal()-minnll);
                minnll = nll->getVal(); 
                rbest = poi.getVal(); 
                rval = rbest;
                poi.setConstant(true);
            }
        }
#endif
  }
  if (points)
    g_outputFile->WriteTObject(points);
  return std::copysign(thisnll - minnll, rbest);
}

double Significance::significanceFromScan(
    RooAbsPdf &pdf, RooAbsData &data, RooRealVar &poi, const RooArgSet *nuisances, double tolerance, int steps) const {
  std::unique_ptr<RooAbsReal> nll(pdf.createNLL(data, RooFit::Constrain(*nuisances)));
  double maxScan = poi.getMax() * 0.7;
  bool stepDown = (bfAlgo_.find("stepDown") != std::string::npos);
  bool twice = (bfAlgo_.find("Twice") != std::string::npos);
  poi.setConstant(false);
  poi.setVal(0.05 * poi.getMax());
  CascadeMinimizer minim0(*nll, CascadeMinimizer::Unconstrained, &poi);
  minim0.setStrategy(0);
  minim0.minimize(g_verbose - 2);
  if (!stepDown) {
    if (poi.getVal() < 0) {
      printf("Minimum found at %s = %8.5f < 0: significance will be zero\n", poi.GetName(), poi.getVal());
      return 0;
    }
    maxScan = poi.getVal() * 1.4;
  } else {
    poi.setVal(0);
  }
  poi.setConstant(true);
  CascadeMinimizer minim(*nll, CascadeMinimizer::Constrained);
  if (!minim.minimize(g_verbose - 2)) {
    std::cerr << "Initial minimization failed. Aborting." << std::endl;
    return -1;
  } else if (g_verbose > 0) {
    printf("Minimum found at %s = %8.5f\n", poi.GetName(), poi.getVal());
  }
  MinimizerSentry minimizerConfig(minimizerAlgoForBF_, minimizerToleranceForBF_);
  std::unique_ptr<RooFitResult> start(minim.save());
  double minnll = nll->getVal(), thisnll = minnll, refnll = thisnll, maxnll = thisnll;
  double rbest = poi.getVal(), rval = rbest;
  TGraph *points = new TGraph(steps + 1);
  points->SetName(Form("nll_scan_%g", g_mass));
  TF1 *fit = new TF1("fit", "[0]*pow(abs(x-[1]), [2])+[3]", 0, poi.getMax());
  fit->SetParNames("norm", "bestfit", "power", "offset");
  points->SetPoint(0, rval, 0);
  if (g_verbose) {
    printf("  %-6s  delta(NLL)\n", poi.GetName());
    printf("%8.5f  %8.5f\n", rval, 0.);
    fflush(stdout);
  }
  for (int i = 1; i < steps; ++i) {
    rval = (maxScan * (stepDown ? i : steps - i - 1)) / steps;
    poi.setVal(rval);
    bool success = minim.improve(g_verbose - 2, /*cascade=*/false);
    thisnll = nll->getVal();
    if (success == false)
      std::cerr << "Minimization failed at " << poi.getVal() << "." << std::endl;
    if (g_verbose) {
      printf("%8.5f  %8.5f\n", rval, thisnll - refnll);
      fflush(stdout);
    }
    points->SetPoint(i, rval, thisnll - refnll);
    if (thisnll < minnll) {
      minnll = thisnll;
      rbest = rval;
    }
    if (rval <= rbest && thisnll > maxnll) {
      maxnll = thisnll;
    }
  }
  if (twice) {
    if (g_verbose) {
      printf("\nPlay it again, sam.\n");
      printf("  %-6s  delta(NLL)\n", poi.GetName());
      fflush(stdout);
    }
    for (int i = steps - 1; i >= 0; --i) {
      rval = (maxScan * (stepDown ? i : steps - i - 1)) / steps;
      if (i == 0 && !stepDown)
        rval = rbest;
      poi.setVal(rval);
      bool success = minim.improve(g_verbose - 2, /*cascade=*/false);
      thisnll = nll->getVal();
      if (success == false)
        std::cerr << "Minimization failed at " << poi.getVal() << "." << std::endl;
      if (g_verbose) {
        printf("%8.5f  %8.5f\n", rval, thisnll - refnll);
        fflush(stdout);
      }
      points->SetPoint(i, rval, thisnll - refnll);
      if (thisnll < minnll) {
        minnll = thisnll;
        rbest = rval;
      }
      if (rval <= rbest && thisnll > maxnll) {
        maxnll = thisnll;
      }
    }
  }
  fit->SetParameters(1, rbest, 2.0, minnll - refnll);
  points->Sort();
  double ret = (maxnll - minnll);
  TFitResultPtr res;
  {
    //MinimizerSentry minimizerConfig(minimizerAlgo_, minimizerTolerance_);
    res = points->Fit(fit, "S0");
  }
  g_outputFile->WriteTObject(points);
  g_outputFile->WriteTObject(fit);
  if (res.Get()->Status() == 0) {
    std::cout << "Using first and last value, the result would be " << ret << std::endl;
    ret = fit->Eval(0) - fit->Eval(fit->GetParameter(1));
    std::cout << "Using output of fit, the result is " << ret << std::endl;
  } else {
    std::cout << "Using first and last value" << std::endl;
  }
  return ret;
}
