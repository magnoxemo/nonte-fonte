#include "Domain.h"
#include "TallyBase.h"
#include "MonteCarloSimulation.h"
#include "FETBase.h"
#include "LegendreFET.h"
#include "HistogramTally.h"


#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>


double beta_pdf(const std::vector<double>& x) {
  double val = x[0];
  if (val < 0.0 || val > 1.0) return 0.0;
  return 6.0 * val * (1.0 - val);
}

int main() {

  nonte_fonte::Domain global_domain(0.0, 1.0);
  int num_bins = 10;
  long n_samples = 100000;
  std::vector<double> mesh;
  for (int i = 0; i <= num_bins; ++i) {
    mesh.push_back(i / static_cast<double>(num_bins));
  }

  nonte_fonte::HistogramTally1D histogram(global_domain, mesh);
  nonte_fonte::LegendreFET tally1(nonte_fonte::Domain(0.0, 0.8), {3});
  nonte_fonte::LegendreFET tally2(nonte_fonte::Domain(0.2, 1.0), {3});

  nonte_fonte::MCSimulation mc(global_domain, beta_pdf, n_samples);

  mc.addTally(&histogram);
  mc.addTally(&tally1);
  mc.addTally(&tally2);

  mc.run();


  return 0;
}