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


double gaussian_2d(const std::vector<double>& p) {

    auto x = p[0];
    auto y = p[1];
    return std::abs(2*x*x*x - y*y + x*x*y -4*x*y*y +5*x*y -3*x +5*y);
}

void write_output(std::ofstream& file, const nonte_fonte::LegendreFET& basis){

  file << "x,y,true_pdf,legendre_pdf\n";

  for (double x = -1.0; x <= 1.0; x += 0.1) {
    for (double y = -1.0; y <= 0; y += 0.1) {
      std::vector<double> point = {x, y};
      file << x << "," << y << "," << gaussian_2d(point) << "," << basis.reconstruct(point) << "\n";
    }
  }
  file.close();
}

int main() {


  nonte_fonte::Domain whole_domain(-2.0, 2.0, -1.0, 1.0);

  nonte_fonte::Domain domain_1(-2.0, 2.0, -1.0, 0);
  nonte_fonte::Domain domain_2(-2.0, 2.0, 0, 1.0);
  nonte_fonte::LegendreFET legendre_1(domain_1, {5,5});
  nonte_fonte::LegendreFET legendre_2(domain_2, {5,5});

  constexpr int n_samples = 500000;

  std::vector<double> x_edges, y_edges;
  for (int i = 0; i <= 10; ++i) {
    double val = -2.0 + i * 4.0 / 10.0;
    x_edges.push_back(val);
    y_edges.push_back(val);
  }
  nonte_fonte::MCSimulation mc(whole_domain, gaussian_2d, n_samples);

  mc.addTally(&legendre_1);
  mc.addTally(&legendre_2);

  mc.run();


  std::ofstream csv_file("example_2d_results.csv");
  write_output(csv_file, legendre_1);


  return 0;
}