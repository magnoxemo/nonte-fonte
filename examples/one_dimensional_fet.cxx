/**
 * @file example_2d.cpp
 * @brief 2D Legendre FET and histogram example
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

// 2D Gaussian
double gaussian_2d(const std::vector<double>& x) {
  double r_sq = x[0]*x[0] + x[1]*x[1];
  return std::exp(-r_sq) / M_PI;
}

int main() {
  std::cout << "========================================\n";
  std::cout << "2D Example\n";
  std::cout << "========================================\n\n";

  // Domain
  Domain domain(-2.0, 2.0, -2.0, 2.0);

  // 2D Legendre FET
  LegendreFET legendre(domain, {8, 8});  // orders (8, 8)

  // 2D structured histogram
  std::vector<double> x_edges, y_edges;
  for (int i = 0; i <= 10; ++i) {
    double val = -2.0 + i * 4.0 / 10.0;
    x_edges.push_back(val);
    y_edges.push_back(val);
  }
  HistogramTally2D histogram(domain, x_edges, y_edges);

  // Monte Carlo
  long n_samples = 500000;
  MCSimulation mc(domain, gaussian_2d, n_samples);

  mc.addTally(&legendre);
  mc.addTally(&histogram);

  mc.run();

  std::cout << "\n";

  // Test points
  std::cout << "PDF Comparison at Test Points:\n";
  std::cout << std::setw(10) << "x" << std::setw(10) << "y"
            << std::setw(15) << "True" << std::setw(15) << "Legendre"
            << std::setw(15) << "Error\n";
  std::cout << std::string(65, '-') << "\n";

  std::vector<std::vector<double>> test_points = {
          {0.0, 0.0}, {0.5, 0.5}, {-0.5, 0.5}, {1.0, 0.0}, {0.0, 1.0}
  };

  for (const auto& p : test_points) {
    double true_val = gaussian_2d(p);
    double recon_val = legendre.reconstruct(p);
    double error = std::abs(true_val - recon_val);

    std::cout << std::setw(10) << std::setprecision(4) << p[0]
              << std::setw(10) << std::setprecision(4) << p[1]
              << std::setw(15) << std::setprecision(6) << true_val
              << std::setw(15) << std::setprecision(6) << recon_val
              << std::setw(15) << std::setprecision(6) << error << "\n";
  }

  std::cout << "\n";

  // Export
  std::ofstream file("example_2d_results.csv");
  file << "x,y,true_pdf,legendre_pdf\n";

  for (double x = -2.0; x <= 2.0; x += 0.1) {
    for (double y = -2.0; y <= 2.0; y += 0.1) {
      std::vector<double> point = {x, y};
      file << x << "," << y << "," << gaussian_2d(point) << ","
           << legendre.reconstruct(point) << "\n";
    }
  }
  file.close();

  std::cout << "Results exported to example_2d_results.csv\n";
  std::cout << "========================================\n";

  return 0;
}