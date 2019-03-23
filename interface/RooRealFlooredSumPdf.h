/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROOREALFLOOREDSUMPDF
#define ROOREALFLOOREDSUMPDF

#include <iostream>
#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooRealProxy.h"
#include "RooAbsCategory.h"
#include "RooCategoryProxy.h"
#include "RooListProxy.h"
#include "RooAICRegistry.h"
#include "RooObjCacheManager.h"
#include "TH3F.h"
#include "TH1.h"
#include "RooDataHist.h"
#include "RooHistFunc.h"

class RooRealFlooredSumPdf : public RooAbsPdf {
public:
  RooRealFlooredSumPdf();
  RooRealFlooredSumPdf(const char* name, const char* title);
  RooRealFlooredSumPdf(const char* name,
                       const char* title,
                       const RooArgList& funcList,
                       const RooArgList& coefList,
                       Bool_t extended = kFALSE);
  RooRealFlooredSumPdf(const RooRealFlooredSumPdf& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new RooRealFlooredSumPdf(*this, newname); }
  virtual ~RooRealFlooredSumPdf();

  Double_t evaluate() const;
  virtual Bool_t checkObservables(const RooArgSet* nset) const;

  virtual Bool_t forceAnalyticalInt(const RooAbsArg&) const { return kTRUE; }
  Int_t getAnalyticalIntegralWN(RooArgSet& allVars,
                                RooArgSet& numVars,
                                const RooArgSet* normSet,
                                const char* rangeName = 0) const;
  Double_t analyticalIntegralWN(Int_t code, const RooArgSet* normSet, const char* rangeName = 0) const;

  const RooArgList& funcList() const { return _funcList; }
  const RooArgList& coefList() const { return _coefList; }

  void setFloor(Double_t val);

  virtual ExtendMode extendMode() const;

  virtual Double_t expectedEvents(const RooArgSet* nset) const;
  virtual Double_t expectedEvents(const RooArgSet& nset) const {
    // Return expected number of events for extended likelihood calculation
    // which is the sum of all coefficients
    return expectedEvents(&nset);
  }

  void printMetaArgs(std::ostream& os) const;

  virtual std::list<Double_t>* binBoundaries(RooAbsRealLValue& /*obs*/, Double_t /*xlo*/, Double_t /*xhi*/) const;
  virtual std::list<Double_t>* plotSamplingHint(RooAbsRealLValue& /*obs*/, Double_t /*xlo*/, Double_t /*xhi*/) const;
  Bool_t isBinnedDistribution(const RooArgSet& obs) const;

protected:
  class CacheElem : public RooAbsCacheElement {
  public:
    CacheElem(){};
    virtual ~CacheElem(){};
    virtual RooArgList containedArgs(Action) {
      RooArgList ret(_funcIntList);
      ret.add(_funcNormList);
      return ret;
    }
    RooArgList _funcIntList;
    RooArgList _funcNormList;
  };
  mutable RooObjCacheManager _normIntMgr;  // The integration cache manager

  Bool_t _haveLastCoef;

  RooListProxy _funcList;  //  List of component FUNCs
  RooListProxy _coefList;  //  List of coefficients
  TIterator* _funcIter;    //! Iterator over FUNC list
  TIterator* _coefIter;    //! Iterator over coefficient list
  Bool_t _extended;        // Allow use as extended p.d.f.
  Bool_t _doFloor;
  Double_t _floorVal;

private:
  ClassDef(RooRealFlooredSumPdf, 2)  // PDF constructed from a sum of (non-pdf) functions
};

#endif
