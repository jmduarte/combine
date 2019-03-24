#include "combine/TestProposal.h"
#include "combine/DebugProposal.h"
#include "combine/VerticalInterpPdf.h"
#include "combine/VerticalInterpHistPdf.h"
#include "combine/AsymPow.h"
#include "combine/AsymQuad.h"
#include "combine/CombDataSetFactory.h"
#include "combine/TH1Keys.h"
#include "combine/RooSimultaneousOpt.h"
#include "combine/SimpleCacheSentry.h"
#include "combine/th1fmorph.h"
#include "combine/HZZ4L_RooCTauPdf_1D.h"
#include "combine/HZZ4L_RooCTauPdf_1D_Expanded.h"
#include "combine/HZZ4L_RooCTauPdf_2D.h"
#include "combine/HZZ4LRooPdfs.h"
#include "combine/HWWLVJRooPdfs.h"
#include "combine/HZZ2L2QRooPdfs.h"
#include "combine/HGGRooPdfs.h"
#include "combine/HZGRooPdfs.h"
#include "combine/SequentialMinimizer.h"
#include "combine/ProcessNormalization.h"
#include "combine/RooRealFlooredSumPdf.h"
#include "combine/RooSpline1D.h"
#include "combine/RooSplineND.h"
#include "combine/RooScaleLOSM.h"
#include "combine/rVrFLikelihood.h"
#include "combine/RooMultiPdf.h"
#include "combine/RooBernsteinFast.h"
#include "combine/SimpleGaussianConstraint.h"
#include "combine/SimplePoissonConstraint.h"
#include "combine/AtlasPdfs.h"
#include "combine/FastTemplateFunc.h"
#include "combine/HZZ4L_RooSpinZeroPdf.h"
#include "combine/HZZ4L_RooSpinZeroPdf_1D.h"
#include "combine/HZZ4L_RooSpinZeroPdf_2D.h"
#include "combine/HZZ4L_RooSpinZeroPdf_phase.h"
#include "combine/VBFHZZ4L_RooSpinZeroPdf.h"
#include "combine/VBFHZZ4L_RooSpinZeroPdf_fast.h"

#include "combine/HZZ4L_RooSpinZeroPdf_1D_fast.h"
#include "combine/HZZ4L_RooSpinZeroPdf_2D_fast.h"
#include "combine/HZZ4L_RooSpinZeroPdf_phase_fast.h"
#include "combine/VVHZZ4L_RooSpinZeroPdf_1D_fast.h"

#include "combine/HWWLVJJRooPdfs.h"
#include "combine/RooMomentMorphND.h"
#include "combine/RooMorphingPdf.h"
#include "combine/RooParametricHist.h"
#include "combine/RooParametricShapeBinPdf.h"
#include "combine/GaussExp.h"
#include "combine/RooDoubleCBFast.h"
#include "combine/CMSHistFunc.h"
#include "combine/CMSHistErrorPropagator.h"
#include "combine/CMSHistFuncWrapper.h"

#include "combine/RooPiecewisePolynomial.h"

#include "combine/RooNCSplineCore.h"
#include "combine/RooNCSpline_1D_fast.h"
#include "combine/RooNCSpline_2D_fast.h"
#include "combine/RooNCSpline_3D_fast.h"
#include "combine/RooFuncPdf.h"

namespace {
  struct dictionary {
    RooBernsteinFast<1> my_RooBernsteinFast_1;
    RooBernsteinFast<2> my_RooBernsteinFast_2;
    RooBernsteinFast<3> my_RooBernsteinFast_3;
    RooBernsteinFast<4> my_RooBernsteinFast_4;
    RooBernsteinFast<5> my_RooBernsteinFast_5;
    RooBernsteinFast<6> my_RooBernsteinFast_6;
    RooBernsteinFast<7> my_RooBernsteinFast_7;
  };
}  // namespace
