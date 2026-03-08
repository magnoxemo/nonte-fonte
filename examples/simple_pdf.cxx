/**
 * @file histogram_example.cpp
 * @brief Histogram and FET combined example matching Python
 */

#include "Domain.h"
#include "TallyBase.h"
#include "MonteCarloSimulation.h"

// FET tallies
#include "FETBase.h"
#include "LegendreFET.h"

// Histogram tallies
#include "HistogramTally.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace local_fet;

// Simple PDF: f(x) = 6x(1-x) on [0,1]
double beta_pdf(const std::vector<double>& x) {
  double val = x[0];
  if (val < 0.0 || val > 1.0) return 0.0;
  return 6.0 * val * (1.0 - val);
}

int main() {
  std::cout << "========================================\n";
  std::cout << "Histogram + FET Example\n";
  std::cout << "========================================\n\n";

  // Global domain
  Domain global_domain(0.0, 1.0);

  // Create histogram with uniform mesh
  int num_bins = 10;
  std::vector<double> mesh;
  for (int i = 0; i <= num_bins; ++i) {
    mesh.push_back(i / static_cast<double>(num_bins));
  }
  HistogramTally1D histogram(global_domain, mesh);

  // Create FET tallies on overlapping domains
  LegendreFET tally1(Domain(0.0, 0.8), {3});
  LegendreFET tally2(Domain(0.2, 1.0), {3});

  // Monte Carlo
  long n_samples = 100000;
  MCSimulation mc(global_domain, beta_pdf, n_samples);

  mc.addTally(&histogram);
  mc.addTally(&tally1);
  mc.addTally(&tally2);

  mc.run();

  std::cout << "\n";

  // Display histogram
  std::cout << "Histogram Results:\n";
  auto centers = histogram.binCenters();
  auto values = histogram.values();
  for (size_t i = 0; i < centers.size(); ++i) {
    std::cout << "  Bin " << i << " (center=" << std::setprecision(3) << centers[i]
              << "): " << std::setprecision(6) << values[i] << "\n";
  }

  std::cout << "\n";

  // Export results
  std::ofstream file("histogram_example_results.csv");
  file << "x,true_pdf,histogram,tally1,tally2\n";

  for (double x = 0.0; x <= 1.0; x += 0.01) {
    std::vector<double> point = {x};
    double true_val = beta_pdf(point);

    // Find histogram value
    double hist_val = 0.0;
    for (size_t i = 0; i < mesh.size() - 1; ++i) {
      if (x >= mesh[i] && x < mesh[i+1]) {
        hist_val = values[i];
        break;
      }
    }

    double tally1_val = tally1.domain().contains(point) ? tally1.reconstruct(point) : 0.0;
    double tally2_val = tally2.domain().contains(point) ? tally2.reconstruct(point) : 0.0;

    file << x << "," << true_val << "," << hist_val << ","
         << tally1_val << "," << tally2_val << "\n";
  }
  file.close();

  std::cout << "Results exported to histogram_example_results.csv\n";
  std::cout << "========================================\n";

  return 0;
}