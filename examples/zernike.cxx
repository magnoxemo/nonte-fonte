/**
 * @file zernike_example.cpp
 * @brief Zernike polynomial example for circular PDF
 */

#include "Domain.h"
#include "TallyBase.h"
#include "MonteCarloSimulation.h"
#include "FETBase.h"
#include "ZernikeFET.h"

#include <iostream>
#include <fstream>
#include <cmath>

using namespace local_fet;

// Gaussian on circular disk
double gaussian_circular(const std::vector<double>& x) {
    double r_sq = x[0]*x[0] + x[1]*x[1];
    double radius = 1.5;

    if (std::sqrt(r_sq) > radius) return 0.0;

    return std::exp(-r_sq) / (M_PI * radius * radius);
}

int main() {
    std::cout << "========================================\n";
    std::cout << "Zernike FET Example\n";
    std::cout << "========================================\n\n";

    // Domain containing circle
    double radius = 1.5;
    Domain domain(-radius, radius, -radius, radius);

    // Zernike tally
    ZernikeFET tally(domain, radius, 8);  // max order = 8

    std::cout << "Number of Zernike modes: " << tally.modes().size() << "\n\n";

    // Monte Carlo
    long n_samples = 500000;
    MCSimulation mc(domain, gaussian_circular, n_samples);
    mc.addTally(&tally);
    mc.run();

    std::cout << "\n";

    // Export results
    std::ofstream file("zernike_results.csv");
    file << "x,y,true_pdf,zernike_pdf\n";

    for (double x = -radius; x <= radius; x += 0.05) {
        for (double y = -radius; y <= radius; y += 0.05) {
            std::vector<double> point = {x, y};
            double r = std::sqrt(x*x + y*y);

            if (r <= radius) {
                file << x << "," << y << ","
                     << gaussian_circular(point) << ","
                     << tally.reconstruct(point) << "\n";
            }
        }
    }
    file.close();

    std::cout << "Results exported to zernike_results.csv\n";
    std::cout << "Plot with: python plotting_utilities.py zernike_results.csv\n";
    std::cout << "========================================\n";

    return 0;
}