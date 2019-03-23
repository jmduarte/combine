#include "TestProposal.h"
#include "DebugProposal.h"
#include "VerticalInterpPdf.h"
#include "VerticalInterpHistPdf.h"
#include "AsymPow.h"
#include "AsymQuad.h"
#include "CombDataSetFactory.h"
#include "TH1Keys.h"
#include "RooSimultaneousOpt.h"
#include "SimpleCacheSentry.h"
#include "th1fmorph.h"
#include "HZZ4L_RooCTauPdf_1D.h"
#include "HZZ4L_RooCTauPdf_1D_Expanded.h"
#include "HZZ4L_RooCTauPdf_2D.h"
#include "HZZ4LRooPdfs.h"
#include "HWWLVJRooPdfs.h"
#include "HZZ2L2QRooPdfs.h"
#include "HGGRooPdfs.h"
#include "HZGRooPdfs.h"
#include "SequentialMinimizer.h"
#include "ProcessNormalization.h"
#include "RooRealFlooredSumPdf.h"
#include "RooSpline1D.h"
#include "RooSplineND.h"
#include "RooScaleLOSM.h"
#include "rVrFLikelihood.h"
#include "RooMultiPdf.h"
#include "RooBernsteinFast.h"
#include "SimpleGaussianConstraint.h"
#include "SimplePoissonConstraint.h"
#include "AtlasPdfs.h"
#include "FastTemplateFunc.h"
#include "HZZ4L_RooSpinZeroPdf.h"
#include "HZZ4L_RooSpinZeroPdf_1D.h"
#include "HZZ4L_RooSpinZeroPdf_2D.h"
#include "HZZ4L_RooSpinZeroPdf_phase.h"
#include "VBFHZZ4L_RooSpinZeroPdf.h"
#include "VBFHZZ4L_RooSpinZeroPdf_fast.h"

#include "HZZ4L_RooSpinZeroPdf_1D_fast.h"
#include "HZZ4L_RooSpinZeroPdf_2D_fast.h"
#include "HZZ4L_RooSpinZeroPdf_phase_fast.h"
#include "VVHZZ4L_RooSpinZeroPdf_1D_fast.h"

#include "HWWLVJJRooPdfs.h"
#include "RooMomentMorphND.h"
#include "RooMorphingPdf.h"
#include "RooParametricHist.h"
#include "RooParametricShapeBinPdf.h"
#include "GaussExp.h"
#include "RooDoubleCBFast.h"
#include "CMSHistFunc.h"
#include "CMSHistErrorPropagator.h"
#include "CMSHistFuncWrapper.h"

#include "RooPiecewisePolynomial.h"

#include "RooNCSplineCore.h"
#include "RooNCSpline_1D_fast.h"
#include "RooNCSpline_2D_fast.h"
#include "RooNCSpline_3D_fast.h"
#include "RooFuncPdf.h"

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
}
