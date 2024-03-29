/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

// Gaussian core + exponential tail
// Souvik Das
// 8/1/2013

#include "Riostream.h"

#include "combine/GaussExp.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include <math.h>
#include "TMath.h"

ClassImp(GaussExp)

    GaussExp::GaussExp(
        const char* name, const char* title, RooAbsReal& _x, RooAbsReal& _p0, RooAbsReal& _p1, RooAbsReal& _p2)
    : RooAbsPdf(name, title),
      x("x", "x", this, _x),
      p0("p0", "p0", this, _p0),
      p1("p1", "p1", this, _p1),
      p2("p2", "p2", this, _p2) {}

GaussExp::GaussExp(const GaussExp& other, const char* name)
    : RooAbsPdf(other, name),
      x("x", this, other.x),
      p0("p0", this, other.p0),
      p1("p1", this, other.p1),
      p2("p2", this, other.p2) {}

Double_t GaussExp::evaluate() const {
  Double_t std = (x - p0) / p1;
  Double_t result = 0;

  if (std < p2) {
    result = exp(-0.5 * pow(std, 2));
  } else {
    result = exp(p2 * p2 / 2. - p2 * std);
  }

  return result;
}
