#include "Domain.h"
#include "TallyBase.h"
#include "MonteCarloSimulation.h"
#include "FETBase.h"
#include "ZernikeFET.h"

#include <iostream>
#include <fstream>
#include <cmath>


double gaussian_circular(const std::vector<double>& x) {
    double r_sq = x[0]*x[0] + x[1]*x[1];
    double radius = 1.5;

    if (std::sqrt(r_sq) > radius) return 0.0;

    return 2+std::exp(-r_sq) / (M_PI * radius * radius);
}

void write_output(std::ofstream& file,double radius, const nonte_fonte::ZernikeFET& basis){

    file << "x,y,true_pdf,legendre_pdf\n";

    for (double x = -radius; x <= radius; x += 0.1) {
        for (double y = -radius; y <= radius; y += 0.1) {
            std::vector<double> point = {x, y};
            double r = std::sqrt(x*x + y*y);
            if (r <= radius) {
                file << x << "," << y << ","
                     << gaussian_circular(point) << ","
                     << basis.reconstruct(point) << "\n";
            }
        }
    }
    file.close();
}

int main() {

    double radius = 1.5;
    long n_samples = 5000000;
    nonte_fonte::Domain domain(-radius, radius, -radius, radius);
    nonte_fonte::ZernikeFET basis(domain, radius, 8);

    std::cout << "Number of Zernike modes: " << basis.modes().size() << "\n\n";

    nonte_fonte::MCSimulation mc(domain, gaussian_circular, n_samples);
    mc.addTally(&basis);
    mc.run();

    std::cout << "\n";

    std::ofstream file("zernike_results.csv");

    write_output(file, radius, basis);


    return 0;
}