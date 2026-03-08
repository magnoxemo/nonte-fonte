/**
 * @file domain_decomposition_example.cpp
 * @brief Domain decomposition example matching Python version
 *
 * Demonstrates multiple FET tallies on different sub-domains.
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


// Grisheimer's pdf
double piecewise_pdf(const std::vector<double>& x) {
    constexpr double norm = 1.0 / 1.51985;
    double val = x[0];

    if (val >= -1.0 && val < -0.5) {
        return norm * std::cos(val) * std::exp(2.0 * val + 1.0);
    } else if (val >= -0.5 && val < 0.5) {
        return norm * std::cos(val);
    } else if (val >= 0.5 && val <= 1.0) {
        return norm * std::cos(val) * std::exp(-val / 2.0 + 0.25);
    }
    return 0.0;
}

int main() {

    nonte_fonte::Domain global_domain(-1.0, 1.0);

    std::vector<nonte_fonte::Domain> domains = {
            nonte_fonte::Domain(-1.0, -0.5),
            nonte_fonte::Domain(-0.5, 0.5),
            nonte_fonte::Domain(0.5, 1.0)
    };
    std::vector<int> orders = {2, 4, 2};

    std::vector<nonte_fonte::LegendreFET> tallies;
    for (size_t i = 0; i < domains.size(); ++i) {
        tallies.emplace_back(domains[i], std::vector<int>{orders[i]});
        std::cout << "Tally " << i+1 << ": ["
                  << domains[i].bounds(0).first << ", "
                  << domains[i].bounds(0).second << "], order=" << orders[i] << "\n";
    }

    std::cout << "\n";

    long n_samples = 30000000;
    nonte_fonte::MCSimulation mc(global_domain, piecewise_pdf, n_samples);

    for (auto& tally : tallies) {
        mc.addTally(&tally);
    }

    mc.run();

    std::cout << "\n";

    // Display coefficients
    std::cout << "Legendre Coefficients:\n";
    for (size_t i = 0; i < tallies.size(); ++i) {
        std::cout << "  Tally " << i+1 << ": [";
        const auto& coeffs = tallies[i].coefficients();
        for (size_t j = 0; j < coeffs.size(); ++j) {
            std::cout << std::setprecision(6) << std::setw(10) << coeffs[j];
            if (j < coeffs.size() - 1) std::cout << ", ";
        }
        std::cout << "]\n";
    }

    std::cout << "\n";

    // Compute errors
    std::cout << "L2 Errors:\n";
    for (size_t i = 0; i < tallies.size(); ++i) {
        double error_sum = 0.0;
        auto [min, max] = domains[i].bounds(0);
        int n_test = 1000;
        double dx = (max - min) / n_test;

        for (int j = 0; j < n_test; ++j) {
            double x_val = min + j * dx;
            std::vector<double> point = {x_val};
            double true_val = piecewise_pdf(point);
            double recon_val = tallies[i].reconstruct(point);
            double diff = true_val - recon_val;
            error_sum += diff * diff * dx;
        }

        std::cout << "  Tally " << i+1 << ": " << std::sqrt(error_sum) << "\n";
    }

    std::cout << "\n";

    // Export results
    std::ofstream file("domain_decomposition_results.csv");
    file << "x,true_pdf";
    for (size_t i = 0; i < tallies.size(); ++i) {
        file << ",tally_" << i;
    }
    file << "\n";

    for (double x = -1.0; x <= 1.0; x += 0.01) {
        std::vector<double> point = {x};
        file << x << "," << piecewise_pdf(point);
        for (auto& tally : tallies) {
            if (tally.domain().contains(point)) {
                file << "," << tally.reconstruct(point);
            } else {
                file << ",0.0";
            }
        }
        file << "\n";
    }
    file.close();

    std::cout << "Results exported to domain_decomposition_results.csv\n";
    std::cout << "========================================\n";

    return 0;
}