#include "combine/GenerateOnly.h"
#include "combine/Combine.h"
#include <iostream>

using namespace RooStats;

GenerateOnly::GenerateOnly() : LimitAlgo("GenerateOnly specific options: none") {}

void GenerateOnly::applyOptions() {}

bool GenerateOnly::run(RooWorkspace *w,
                       RooStats::ModelConfig *mc_s,
                       RooStats::ModelConfig *mc_b,
                       RooAbsData &data,
                       double &limit,
                       double &limitErr,
                       const double *hint) {
  if (g_verbose > 0) {
    std::cout << "generate toy samples only; no limit computation " << std::endl;
  }
  // Kill the Fill since there is nothing to do
  Combine::toggleGlobalFillTree(false);
  return true;
}
