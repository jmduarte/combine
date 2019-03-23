#include "combine/interface/TestProposal.h"
#include "combine/interface/DebugProposal.h"
#include "combine/interface/VerticalInterpPdf.h"
#include "combine/interface/VerticalInterpHistPdf.h"
#include "combine/interface/AsymPow.h"
#include "combine/interface/AsymQuad.h"
#include "combine/interface/CombDataSetFactory.h"
#include "combine/interface/TH1Keys.h"
#include "combine/interface/RooSimultaneousOpt.h"
#include "combine/interface/SimpleCacheSentry.h"
#include "combine/interface/th1fmorph.h"
#include "combine/interface/HZZ4L_RooCTauPdf_1D.h"
#include "combine/interface/HZZ4L_RooCTauPdf_1D_Expanded.h"
#include "combine/interface/HZZ4L_RooCTauPdf_2D.h"
#include "combine/interface/HZZ4LRooPdfs.h"
#include "combine/interface/HWWLVJRooPdfs.h"
#include "combine/interface/HZZ2L2QRooPdfs.h"
#include "combine/interface/HGGRooPdfs.h"
#include "combine/interface/HZGRooPdfs.h"
#include "combine/interface/SequentialMinimizer.h"
#include "combine/interface/ProcessNormalization.h"
#include "combine/interface/RooRealFlooredSumPdf.h"
#include "combine/interface/RooSpline1D.h"
#include "combine/interface/RooSplineND.h"
#include "combine/interface/RooScaleLOSM.h"
#include "combine/interface/rVrFLikelihood.h"
#include "combine/interface/RooMultiPdf.h"
#include "combine/interface/RooBernsteinFast.h"
#include "combine/interface/SimpleGaussianConstraint.h"
#include "combine/interface/SimplePoissonConstraint.h"
#include "combine/interface/AtlasPdfs.h"
#include "combine/interface/FastTemplateFunc.h"
#include "combine/interface/HZZ4L_RooSpinZeroPdf.h"
#include "combine/interface/HZZ4L_RooSpinZeroPdf_1D.h"
#include "combine/interface/HZZ4L_RooSpinZeroPdf_2D.h"
#include "combine/interface/HZZ4L_RooSpinZeroPdf_phase.h"
#include "combine/interface/VBFHZZ4L_RooSpinZeroPdf.h"
#include "combine/interface/VBFHZZ4L_RooSpinZeroPdf_fast.h"

#include "combine/interface/HZZ4L_RooSpinZeroPdf_1D_fast.h"
#include "combine/interface/HZZ4L_RooSpinZeroPdf_2D_fast.h"
#include "combine/interface/HZZ4L_RooSpinZeroPdf_phase_fast.h"
#include "combine/interface/VVHZZ4L_RooSpinZeroPdf_1D_fast.h"

#include "combine/interface/HWWLVJJRooPdfs.h"
#include "combine/interface/RooMomentMorphND.h"
#include "combine/interface/RooMorphingPdf.h"
#include "combine/interface/RooParametricHist.h"
#include "combine/interface/RooParametricShapeBinPdf.h"
#include "combine/interface/GaussExp.h"
#include "combine/interface/RooDoubleCBFast.h"
#include "combine/interface/CMSHistFunc.h"
#include "combine/interface/CMSHistErrorPropagator.h"
#include "combine/interface/CMSHistFuncWrapper.h"

#include "combine/interface/RooPiecewisePolynomial.h"

#include "combine/interface/RooNCSplineCore.h"
#include "combine/interface/RooNCSpline_1D_fast.h"
#include "combine/interface/RooNCSpline_2D_fast.h"
#include "combine/interface/RooNCSpline_3D_fast.h"
#include "combine/interface/RooFuncPdf.h"

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
