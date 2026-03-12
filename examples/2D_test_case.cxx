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


double gaussian_2d(const std::vector<double>& x) {
  double r_sq = x[0]*x[0] + x[1]*x[1];
  return std::exp(-r_sq) / M_PI;
}

void write_output(std::ofstream& file, const nonte_fonte::LegendreFET& basis){

  file << "x,y,true_pdf,legendre_pdf\n";

  for (double x = -2.0; x <= 2.0; x += 0.1) {
    for (double y = -2.0; y <= 2.0; y += 0.1) {
      std::vector<double> point = {x, y};
      file << x << "," << y << "," << gaussian_2d(point) << ","
           << basis.reconstruct(point) << "\n";
    }
  }
  file.close();
}

int main() {

  nonte_fonte::Domain domain(-2.0, 2.0, -2.0, 2.0);
  nonte_fonte::LegendreFET legendre(domain, {8, 8});  // orders (8, 8)
  constexpr int n_samples = 500000;

  std::vector<double> x_edges, y_edges;
  for (int i = 0; i <= 10; ++i) {
    double val = -2.0 + i * 4.0 / 10.0;
    x_edges.push_back(val);
    y_edges.push_back(val);
  }
  nonte_fonte::HistogramTally2D histogram(domain, x_edges, y_edges);
  nonte_fonte::MCSimulation mc(domain, gaussian_2d, n_samples);

  mc.addTally(&legendre);
  mc.addTally(&histogram);

  mc.run();

  std::cout << "\n";

  std::cout << "PDF Comparison at Test Points:\n";
  std::cout << std::setw(10) << "x" << std::setw(10) << "y"
            << std::setw(15) << "True" << std::setw(15) << "Legendre"
            << std::setw(15) << "Error\n";
  std::cout << std::string(65, '-') << "\n";


  std::ofstream csv_file("example_2d_results.csv");
  write_output(csv_file, legendre);


  return 0;
}