#include "Domain.h"
#include "TallyBase.h"
#include "MonteCarloSimulation.h"
#include "FETBase.h"
#include "LegendreFET.h"
#include "HistogramTally.h"
#include "Utilities.h"


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
    auto domains = nonte_fonte::lin_space_domain(-1, 1, 10);
    std::vector<int> orders = {2,2,10,2,2,2,2,2,2,2};

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


    // write in CSV
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


    return 0;
}