cmake_minimum_required(VERSION 3.2)

project(_combine)

set(SOURCE_FILES
    AsimovUtils.cc
    CloseCoutSentry.cc
    GBRMath.cc
    HZZ4L_RooSpinZeroPdf_1D_fast.cc
    ProfilingTools.cc
    RooNCSplineFactory_1D.cc
    SimpleCacheSentry.cc
    vectorized.cc
    AsymPow.cc
    GenerateOnly.cc
    HZZ4L_RooSpinZeroPdf_2D.cc
    RobustHesse.cc
    RooNCSplineFactory_2D.cc
    SimpleGaussianConstraint.cc
    VectorizedGaussian.cc
    AsymptoticLimits.cc
    CMSHistErrorPropagator.cc
    GoodnessOfFit.cc
    HZZ4L_RooSpinZeroPdf_2D_fast.cc
    RooBernsteinFast.cc
    RooNCSplineFactory_3D.cc
    SimplePoissonConstraint.cc
    vectorized.h
    AsymQuad.cc
    CMSHistFunc.cc
    HGGRooPdfs.cc
    HZZ4L_RooSpinZeroPdf.cc
    RooCheapProduct.cc
    RooParametricHist.cxx
    SimplerLikelihoodRatioTestStat.cc
    VectorizedHistFactoryPdfs.cc
    AtlasPdfs.cxx
    CMSHistFuncWrapper.cc
    HWWLVJJRooPdfs.cc
    HZZ4L_RooSpinZeroPdf_phase.cc
    RooDoubleCBFast.cc
    RooParametricShapeBinPdf.cc
    SimplerLikelihoodRatioTestStatExt.cc
    VectorizedSimplePdfs.cc
    BayesianFlatPrior.cc
    CombDataSetFactory.cc
    HWWLVJRooPdfs.cxx
    HZZ4L_RooSpinZeroPdf_phase_fast.cc
    RooFitGlobalKillSentry.cc
    RooPiecewisePolynomial.cc
    tdrstyle.cc
    VerticalInterpHistPdf.cc
    BayesianToyMC.cc
    Combine.cc
    HybridNew.cc
    JacknifeQuantile.cc
    RooFuncPdf.cc
    RooRealFlooredSumPdf.cc
    TestProposal.cc
    VerticalInterpPdf.cc
    BestFitSigmaTestStat.cc
    DebugProposal.cc
    HZGRooPdfs.cxx
    Logger.cc
    RooMomentMorphND.cxx
    RooScaleLOSM.cc
    th1fmorph.cc
    VVHZZ4L_RooSpinZeroPdf_1D_fast.cc
    CachingMultiPdf.cc
    FastTemplate_Old.cc
    HZZ2L2QRooPdfs.cc
    MarkovChainMC.cc
    RooMorphingPdf.cc
    RooSimultaneousOpt.cc
    TH1Keys.cc
    CachingNLL.cc
    FeldmanCousins.cc
    HZZ4L_RooCTauPdf_1D.cc
    MultiDimFit.cc
    RooMultiPdf.cxx
    RooSpline1D.cc
    ToyMCSamplerOpt.cc
    CascadeMinimizer.cc
    FitDiagnostics.cc
    HZZ4L_RooCTauPdf_1D_Expanded.cc
    PdfDiagonalizer.cc
    RooNCSpline_1D_fast.cc
    RooSplineND.cc
    utils.cc
    ChannelCompatibilityCheck.cc
    FitterAlgoBase.cc
    HZZ4L_RooCTauPdf_2D.cc
    ProcessNormalization.cc
    RooNCSpline_2D_fast.cc
    rVrFLikelihood.cc
    VBFHZZ4L_RooSpinZeroPdf.cc
    FnTimer.cc
    HZZ4LRooPdfs.cc
    ProfiledLikelihoodRatioTestStat.cc
    RooNCSpline_3D_fast.cc
    SequentialMinimizer.cc
    VBFHZZ4L_RooSpinZeroPdf_fast.cc
    GaussExp.cxx
    HZZ4L_RooSpinZeroPdf_1D.cc
    ProfiledLikelihoodRatioTestStatExt.cc
    RooNCSplineCore.cc
    Significance.cc
    VectorizedCB.cc
)

find_package(ROOT REQUIRED COMPONENTS Core Imt RIO Net Hist Graf Graf3d Gpad ROOTVecOps Tree TreePlayer Rint Postscript
                                      Matrix Physics MathCore Thread MultiProc ROOTDataFrame MathMore RooFit RooFitCore
                                      RooStats Minuit Foam HistFactory)
include(${ROOT_USE_FILE})

ROOT_GENERATE_DICTIONARY(combine_xr ${COMBINE_HEADERS_DIR}/combine/classes.h
                         LINKDEF ${COMBINE_HEADERS_DIR}/combine/classesLinkDef.h)
add_library(Combine SHARED ${SOURCE_FILES} combine_xr.cxx)
target_link_libraries(Combine PUBLIC ${ROOT_LIBRARIES})

add_subdirectory(pybind11)

pybind11_add_module(_combine combine.cpp)
target_link_libraries(_combine PUBLIC Combine)

install(FILES ${PROJECT_BINARY_DIR}/src/libcombine_xr_rdict.pcm DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(TARGETS _combine DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
